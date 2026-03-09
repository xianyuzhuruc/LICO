#pragma once

#if defined(__linux__)

#include <chrono>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <string>
#include <vector>

#include <asm/unistd.h>
#include <linux/perf_event.h>
#include <sys/ioctl.h>
#include <unistd.h>

namespace _perf_sysfs_ {

inline bool read_all(const std::string& path, std::string& out) {
  std::ifstream ifs(path);
  if (!ifs) return false;
  out.assign(std::istreambuf_iterator<char>(ifs), std::istreambuf_iterator<char>());
  return true;
}
inline std::string trim(const std::string& s) {
  size_t b = s.find_first_not_of(" \t\r\n");
  size_t e = s.find_last_not_of(" \t\r\n");
  if (b == std::string::npos) return "";
  return s.substr(b, e - b + 1);
}
inline std::map<std::string, uint64_t> parse_kv(const std::string& spec) {
  std::map<std::string, uint64_t> kv;
  std::string s = spec;
  for (auto& ch : s) if (ch == '\n') ch = ',';
  size_t pos = 0;
  while (pos < s.size()) {
    size_t next = s.find(',', pos);
    std::string item = trim(s.substr(pos, (next == std::string::npos ? s.size() : next) - pos));
    if (!item.empty()) {
      size_t eq = item.find('=');
      std::string key = trim(item.substr(0, eq == std::string::npos ? item.size() : eq));
      uint64_t val = 1;
      if (eq != std::string::npos) {
        std::string v = trim(item.substr(eq + 1));
        if (v.size() > 2 && v[0] == '0' && (v[1] == 'x' || v[1] == 'X')) {
          val = std::stoull(v, nullptr, 16);
        } else {
          val = std::stoull(v, nullptr, 0);
        }
      }
      kv[key] = val;
    }
    if (next == std::string::npos) break;
    pos = next + 1;
  }
  return kv;
}

struct FieldMap {
  std::string target; // "config" / "config1" / "config2"
  int lo = 0, hi = 0;
};
inline bool parse_format_line(const std::string& line, FieldMap& m) {
  auto s = trim(line);
  auto pos = s.find(':');
  if (pos == std::string::npos) return false;
  m.target = s.substr(0, pos);
  auto range = s.substr(pos + 1);
  auto dash = range.find('-');
  if (dash == std::string::npos) {
    m.lo = std::stoi(range);
    m.hi = m.lo;
  } else {
    m.lo = std::stoi(range.substr(0, dash));
    m.hi = std::stoi(range.substr(dash + 1));
  }
  return true;
}
inline void apply_field(uint64_t& target, const FieldMap& f, uint64_t value) {
  int width = f.hi - f.lo + 1;
  uint64_t mask = (width >= 64) ? ~0ULL : ((1ULL << width) - 1ULL);
  uint64_t v = (value & mask) << f.lo;
  uint64_t clearMask = ~(mask << f.lo);
  target = (target & clearMask) | v;
}

struct AliasResolved {
  bool ok = false;
  std::string pmu_root;   // e.g. /sys/devices/cpu_core
  std::string pmu_tag;    // "core" / "atom" / "cpu"
  uint32_t pmu_type = 0;
  uint64_t config = 0, config1 = 0, config2 = 0;
};

// 在所有可能 PMU 根目录下解析 alias，返回所有成功解析的条目（可能有多个 PMU）
inline std::vector<AliasResolved> resolve_alias_all(const std::string& alias) {
  const std::vector<std::pair<std::string,std::string>> roots = {
    {"/sys/devices/cpu",      "cpu"},
    {"/sys/bus/event_source/devices/cpu",      "cpu"},
    {"/sys/devices/cpu_core", "core"},
    {"/sys/bus/event_source/devices/cpu_core", "core"},
    {"/sys/devices/cpu_atom", "atom"},
    {"/sys/bus/event_source/devices/cpu_atom", "atom"}
  };
  const std::vector<std::string> variants = {
    alias,
    [] (std::string s){ for (auto& c: s) if (c=='.') c='-'; return s; }(alias),
    [] (std::string s){ for (auto& c: s) if (c=='.') c='_'; return s; }(alias)
  };

  std::vector<AliasResolved> out;
  for (const auto& [root, tag] : roots) {
    std::ifstream f1(root + "/format"), f2(root + "/events"), f3(root + "/type");
    if (!f1.good() || !f2.good() || !f3.good()) continue;

    std::string spec;
    bool found = false;
    for (const auto& v : variants) {
      if (read_all(root + "/events/" + v, spec)) { found = true; break; }
    }
    if (!found) continue;

    AliasResolved ar;
    std::string type_str;
    if (!read_all(root + "/type", type_str)) continue;
    ar.pmu_type = static_cast<uint32_t>(std::stoul(trim(type_str)));
    ar.pmu_root = root;
    ar.pmu_tag  = tag;

    auto kv = parse_kv(spec);
    for (const auto& [key, val] : kv) {
      std::string fmt;
      if (!read_all(root + "/format/" + key, fmt)) continue;
      std::stringstream ss(fmt);
      std::string line;
      while (std::getline(ss, line)) {
        FieldMap fm;
        if (!parse_format_line(line, fm)) continue;
        if (fm.target == "config")  apply_field(ar.config,  fm, val);
        if (fm.target == "config1") apply_field(ar.config1, fm, val);
        if (fm.target == "config2") apply_field(ar.config2, fm, val);
      }
    }
    ar.ok = true;
    out.push_back(ar);
  }
  return out;
}

} // namespace _perf_sysfs_

struct PerfEvent {

   struct event {
      struct read_format {
         uint64_t value;
         uint64_t time_enabled;
         uint64_t time_running;
         uint64_t id;
      };

      perf_event_attr pe;
      int fd;
      bool is_raw = false;
      read_format prev;
      read_format data;

      double readCounter() {
         if (is_raw) {
            return static_cast<double>(data.value - prev.value);
         }
         uint64_t delta_running  = data.time_running - prev.time_running;
         uint64_t delta_enabled  = data.time_enabled - prev.time_enabled;
         uint64_t delta_value    = data.value - prev.value;
         if (delta_running == 0) return static_cast<double>(delta_value);
         double mux = static_cast<double>(delta_enabled) / static_cast<double>(delta_running);
         return static_cast<double>(delta_value) * mux;
      }
   };

   std::vector<event> events;
   std::vector<std::string> names;
   std::chrono::time_point<std::chrono::steady_clock> startTime;
   std::chrono::time_point<std::chrono::steady_clock> stopTime;

   PerfEvent() {
      registerCounter("cycles", PERF_TYPE_HARDWARE, PERF_COUNT_HW_CPU_CYCLES);
      registerCounter("instructions", PERF_TYPE_HARDWARE, PERF_COUNT_HW_INSTRUCTIONS);
      registerCounter("L1-misses", PERF_TYPE_HW_CACHE, PERF_COUNT_HW_CACHE_L1D|(PERF_COUNT_HW_CACHE_OP_READ<<8)|(PERF_COUNT_HW_CACHE_RESULT_MISS<<16));
      registerCounter("LLC-misses", PERF_TYPE_HARDWARE, PERF_COUNT_HW_CACHE_MISSES);
      registerCounter("branch-misses", PERF_TYPE_HARDWARE, PERF_COUNT_HW_BRANCH_MISSES);
      registerCounter("task-clock", PERF_TYPE_SOFTWARE, PERF_COUNT_SW_TASK_CLOCK);

      registerCounter("L1D-stores",
          PERF_TYPE_HW_CACHE,
          PERF_COUNT_HW_CACHE_L1D |
              (PERF_COUNT_HW_CACHE_OP_WRITE << 8) |
              (PERF_COUNT_HW_CACHE_RESULT_ACCESS << 16));

      registerCounter("dTLB-misses", PERF_TYPE_HW_CACHE, PERF_COUNT_HW_CACHE_DTLB | (PERF_COUNT_HW_CACHE_OP_READ << 8) | (PERF_COUNT_HW_CACHE_RESULT_MISS << 16));

      registerCounter("LLC-l",
          PERF_TYPE_HW_CACHE,
          PERF_COUNT_HW_CACHE_LL |
              (PERF_COUNT_HW_CACHE_OP_READ << 8) |
              (PERF_COUNT_HW_CACHE_RESULT_ACCESS << 16));

      registerCounter("LLC-l-misses",
          PERF_TYPE_HW_CACHE,
          PERF_COUNT_HW_CACHE_LL |
              (PERF_COUNT_HW_CACHE_OP_READ << 8) |
              (PERF_COUNT_HW_CACHE_RESULT_MISS << 16));

      registerCounter("LLC-w-accesses",
          PERF_TYPE_HW_CACHE,
          PERF_COUNT_HW_CACHE_LL |
              (PERF_COUNT_HW_CACHE_OP_WRITE << 8) |
              (PERF_COUNT_HW_CACHE_RESULT_ACCESS << 16));


      // 打开所有事件
      for (unsigned i=0; i<events.size(); i++) {
         auto& event = events[i];
         event.fd = syscall(__NR_perf_event_open, &event.pe, 0, -1, -1, 0);
         if (event.fd < 0) {
            std::cerr << "Error opening counter " << names[i] << std::endl;
            events.resize(0);
            names.resize(0);
            return;
         }
      }
   }

   void registerRawCounter(const std::string &name, uint64_t config) {
      names.push_back(name);
      events.emplace_back();
      auto &ev = events.back();
      memset(&ev.pe, 0, sizeof(ev.pe));
      memset(&ev.prev, 0, sizeof(ev.prev));
      memset(&ev.data, 0, sizeof(ev.data));

      ev.pe.type   = PERF_TYPE_RAW;     // 将在调用后被设置为具体 PMU type
      ev.pe.size   = sizeof(ev.pe);
      ev.pe.config = config;
      ev.pe.disabled       = 1;
      ev.pe.exclude_kernel = 0;
      ev.pe.exclude_hv     = 0;
      ev.pe.read_format    = 0;         // RAW：只读 u64 value
      ev.is_raw            = true;
   }

   void registerCounter(const std::string& name, uint64_t type, uint64_t eventID) {
      names.push_back(name);
      events.push_back(event());
      auto& event = events.back();
      auto& pe = event.pe;
      memset(&pe, 0, sizeof(struct perf_event_attr));
      memset(&event.prev, 0, sizeof(event.prev));
      memset(&event.data, 0, sizeof(event.data));
      event.is_raw = false;

      pe.type = type;
      pe.size = sizeof(struct perf_event_attr);
      pe.config = eventID;
      pe.disabled = true;
      pe.inherit = 1;
      pe.inherit_stat = 0;
      pe.exclude_kernel = false;
      pe.exclude_hv = false;
      pe.read_format = PERF_FORMAT_TOTAL_TIME_ENABLED | PERF_FORMAT_TOTAL_TIME_RUNNING;
   }

   void startCounters() {
      for (unsigned i=0; i<events.size(); i++) {
         auto& event = events[i];
         ioctl(event.fd, PERF_EVENT_IOC_RESET, 0);
         ioctl(event.fd, PERF_EVENT_IOC_ENABLE, 0);

         if (event.is_raw) {
            uint64_t val = 0;
            ssize_t n = read(event.fd, &val, sizeof(val));
            if (n != sizeof(val)) {
               std::cerr << "Error reading RAW counter " << names[i] << std::endl;
            }
            event.prev.value        = val;
            event.prev.time_enabled = 1;
            event.prev.time_running = 1;
         } else {
            ssize_t n = read(event.fd, &event.prev, sizeof(uint64_t) * 3);
            if (n != ssize_t(sizeof(uint64_t) * 3))
               std::cerr << "Error reading counter " << names[i] << std::endl;
         }
      }
      startTime = std::chrono::steady_clock::now();
   }

   ~PerfEvent() {
      for (auto& event : events) {
         close(event.fd);
      }
   }

   void stopCounters() {
      stopTime = std::chrono::steady_clock::now();
      for (unsigned i=0; i<events.size(); i++) {
         auto& event = events[i];
         if (event.is_raw) {
            uint64_t val = 0;
            ssize_t n = read(event.fd, &val, sizeof(val));
            if (n != sizeof(val)) {
               std::cerr << "Error reading RAW counter " << names[i] << std::endl;
            }
            event.data.value        = val;
            event.data.time_enabled = 1;
            event.data.time_running = 1;
         } else {
            ssize_t n = read(event.fd, &event.data, sizeof(uint64_t) * 3);
            if (n != ssize_t(sizeof(uint64_t) * 3))
               std::cerr << "Error reading counter " << names[i] << std::endl;
         }
         ioctl(event.fd, PERF_EVENT_IOC_DISABLE, 0);
      }
   }

   double getDuration() {
      return std::chrono::duration<double>(stopTime - startTime).count() * 1000 * 1000 * 1000;
   }

   double getIPC() {
      return getCounter("instructions") / getCounter("cycles");
   }

   double getCPUs() {
      return getCounter("task-clock") / (getDuration() * 1e9);
   }

   double getGHz() {
      return getCounter("cycles") / getCounter("task-clock");
   }

   double getMissRate() {
      return getCounter("LLC-l-misses") / getCounter("LLC-l");
   }

   double getCounter(const std::string& name) {
      for (unsigned i=0; i<events.size(); i++)
         if (names[i]==name)
            return events[i].readCounter();
      return -1;
   }

   // 便于合计：求前缀相同的一组计数之和
   double sumCountersByPrefix(const std::string& prefix) {
      double s = 0.0;
      for (unsigned i=0; i<events.size(); i++) {
         if (names[i].rfind(prefix, 0) == 0) { // 以 prefix 开头
           s += events[i].readCounter();
         }
      }
      return s;
   }

   static void printCounter(std::ostream& headerOut, std::ostream& dataOut, std::string name, std::string counterValue,bool addComma=true) {
     auto width=std::max(name.length(),counterValue.length());
     headerOut << std::setw(width) << name << (addComma ? "," : "") << " ";
     dataOut << std::setw(width) << counterValue << (addComma ? "," : "") << " ";
   }

   template <typename T>
   static void printCounter(std::ostream& headerOut, std::ostream& dataOut, std::string name, T counterValue,bool addComma=true) {
     std::stringstream stream;
     stream << std::fixed << std::setprecision(2) << counterValue;
     PerfEvent::printCounter(headerOut,dataOut,name,stream.str(),addComma);
   }

   void printReport(std::ostream& headerOut, std::ostream& dataOut, uint64_t normalizationConstant) {
      if (!events.size())
         return;

      // 原样输出所有已注册计数器
      for (unsigned i=0; i<events.size(); i++) {
         printCounter(headerOut,dataOut,names[i],events[i].readCounter()/normalizationConstant);
      }


      // 你原来的派生指标
      printCounter(headerOut,dataOut,"IPC",getIPC());
      printCounter(headerOut,dataOut,"GHz",getGHz(),false);
   }
};

struct BenchmarkParameters {

  void setParam(const std::string& name,const std::string& value) {
    params[name]=value;
  }

  void setParam(const std::string& name,const char* value) {
    params[name]=value;
  }

  template <typename T>
  void setParam(const std::string& name,T value) {
    setParam(name,std::to_string(value));
  }

  void printParams(std::ostream& header,std::ostream& data) {
    for (auto& p : params) {
      PerfEvent::printCounter(header,data,p.first,p.second);
    }
  }

  BenchmarkParameters(std::string name="") {
    if (name.length())
      setParam("name",name);
  }

  private:
  std::map<std::string,std::string> params;
};

struct PerfEventBlock {
   PerfEvent e;
   uint64_t scale;
   BenchmarkParameters parameters;
   bool printHeader;

   PerfEventBlock(uint64_t scale = 1, BenchmarkParameters params = {}, bool printHeader = true)
       : scale(scale),
         parameters(params),
         printHeader(printHeader) {
     e.startCounters();
   }

   ~PerfEventBlock() {
     e.stopCounters();
     std::stringstream header;
     std::stringstream data;
     parameters.printParams(header,data);
     PerfEvent::printCounter(header,data,"time [us]",e.getDuration());
     e.printReport(header, data, scale);
     if (printHeader)
       std::cout << header.str() << std::endl;
     std::cout << "RESULT: " << data.str() << std::endl;
   }
};

#else
#include <ostream>
struct PerfEvent {
  void startCounters() {}
  void stopCounters() {}
  void printReport(std::ostream&, uint64_t) {}
  template <class T> void setParam(const std::string&, const T&) {};
};

struct BenchmarkParameters {
};

struct PerfEventBlock {
  PerfEventBlock(uint64_t = 1, BenchmarkParameters = {}, bool = true) {};
  PerfEventBlock(PerfEvent e, uint64_t = 1, BenchmarkParameters = {}, bool = true) {};
};
#endif