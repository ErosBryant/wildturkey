//
// Created by daiyi on 2020/02/02.
//


#include "learned_index.h"

#include "db/version_set.h"
#include <cassert>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <sstream>
#include <random>

#include <utility>

#include "util/mutexlock.h"

#include "util.h"
#include "Q_table.h"
#include <cmath> 
#include <chrono>

namespace adgMod {

std::pair<uint64_t, uint64_t> LearnedIndexData::GetPosition(
    const Slice& target_x) const {
  assert(string_segments.size() > 1);
  ++served;

  // check if the key is within the model bounds
  uint64_t target_int = SliceToInteger(target_x);
  if (target_int > max_key) return std::make_pair(size, size);
  if (target_int < min_key) return std::make_pair(size, size);

  // binary search between segments
  uint32_t left = 0, right = (uint32_t)string_segments.size() - 1;
  while (left != right - 1) {
    uint32_t mid = (right + left) / 2;
    if (target_int < string_segments[mid].x)
      right = mid;
    else
      left = mid;
  }

  // calculate the interval according to the selected segment
  // predicted position
  double result =
      target_int * string_segments[left].k + string_segments[left].b;
  result = is_level ? result / 2 : result;
  //double error_bound = this->meta->error;     // 调用对应SST的error_bound
  //std::cout << this->meta->error << std::endl;
  //double error_bound = this->meta->error > 3 ? this->meta->error : 40;     // 调用对应SST的error_bound
  double error_bound = this->meta->error;
  //error_bound = 50;
  uint64_t lower =
      result - error > 0 ? (uint64_t)std::floor(result - error_bound) : 0;
  uint64_t upper = (uint64_t)std::ceil(result + error_bound);
  if (lower >= size) return std::make_pair(size, size);
  upper = upper < size ? upper : size - 1;
  //std::cout << error_bound << std::endl;
  //                printf("%s %s %s\n", string_keys[lower].c_str(),
  //                string(target_x.data(), target_x.size()).c_str(),
  //                string_keys[upper].c_str()); assert(target_x >=
  //                string_keys[lower] && target_x <= string_keys[upper]);
  return std::make_pair(lower, upper);
}

uint64_t LearnedIndexData::MaxPosition() const { return size - 1; }

double LearnedIndexData::GetError() const { return this->meta->error; }

double LearnedIndexData::getRandomAction(double error_) const {
    double epsilon = 0.3; // 探索率 30%

    thread_local std::random_device rd;
    thread_local std::mt19937 gen(rd());
    thread_local std::uniform_real_distribution<> dis(0.0, 1.0);

    double random_value = dis(gen);

    bool is_exploration = random_value < epsilon;
    //std::cout << "Random Value: " << random_value << ", Epsilon: " << epsilon
              //<< ", Exploration: " << (is_exploration ? "Yes" : "No") << std::endl;

    if (is_exploration) {
        // 探索：选择与当前 error_ 相邻的值
        std::vector<double> possible_errors;

        // 定义最小和最大允许的 error_bound
        double min_error = 4;
        double max_error = 40;

        // 添加比当前 error_ 小的值
        if (error_ - 4 >= min_error) {
            possible_errors.push_back(error_ - 4);
        }
        // 添加比当前 error_ 大的值
        if (error_ + 4 <= max_error) {
            possible_errors.push_back(error_ + 4);
        }

        // 确保有候选值
        if (!possible_errors.empty()) {
            std::uniform_int_distribution<> error_dis(0, possible_errors.size() - 1);
            double new_error = possible_errors[error_dis(gen)];
            return new_error;
        } else {
            // 如果没有相邻的候选值，保持不变
            return error_;
        }
    } else {
        // 利用：选择当前最优动作（即当前的 error_）
        return error_;
    }
}




// Actual function doing learning
bool LearnedIndexData::Learn() {

  srand(static_cast<unsigned>(time(0)));  // 设置随机数种子

  // FILL IN GAMMA (error)
  // PLR plr = PLR(error);                           // error可调整为q-table对应的error

  // check if data if filled
  if (string_keys.empty()) assert(false);

  // fill in some bounds for the model
  uint64_t temp = atoll(string_keys.back().c_str());
  min_key = atoll(string_keys.front().c_str());   // SST内最小Key, 为Uint64型
  max_key = atoll(string_keys.back().c_str());    // SST内最大Key, 为Uint64型
  size = string_keys.size();                      // SST内Key的数量, 为Uint64型   以上三变量均可用于计算SST内key分布的模拟密度
  inverse_density = static_cast<uint64_t>((max_key - min_key) / size);            // SST内key分布的模拟密度
  // 预定将inverse_density < 40 为一档， < 70为二档， < 100为3档， >100为4档， 档位越大，密度越小
  // 为了方便计算，将inverse_density的值设定为档位的值，即inverse_density = 1, 2, 3, 4
  if(inverse_density < 10) inverse_density = 0;
  else if(inverse_density < 30) inverse_density = 1;
  else if(inverse_density < 70) inverse_density = 2;
  else inverse_density = 3;
  double temp_error = getRandomAction(adgMod::getQTableManagerInstance().Q_table[inverse_density].error_bound);  // 根据档位选择对应的error_bound
  //error = 32;
  //temp_error = 32;
  PLR plr = PLR(temp_error);                           // error可调整为q-table对应的error
  this->meta = new FileMetaData();
  this->meta->error = temp_error;   
  

  // actual training
  // 记录训练时间
  auto start = std::chrono::high_resolution_clock::now();
  std::vector<Segment> segs = plr.train(string_keys, !is_level);
  auto end = std::chrono::high_resolution_clock::now();
  uint64_t build_time = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();

  // 计算平均模型加载时间、平均校正时间， 参照zipf's law
  int segment_count = segs.size();
  double sum_prob = 0.0;
  for (int i = 1; i <= segment_count; ++i) {
    sum_prob += 1.0 / pow(i, 1);  // 计算 Zipf 总和
  }
  double weighted_sum = 0.0;
  for (int i = 1; i <= segment_count; ++i) {
    double P_i = (1.0 / pow(i, 1)) / sum_prob;  // 段 i 的访问概率
    double S_i = log2(i);  // 段 i 的查找步数
    weighted_sum += P_i * S_i;  // 加权累加
  }
  double reward = adgMod::getQTableManagerInstance().compute_reward(
    inverse_density,
    weighted_sum,    // new_load_model_cost
    temp_error,      // new_correct_time
    build_time       // new_build_time
  );
  double Q_value = adgMod::getQTableManagerInstance().compute_q_value(inverse_density, reward);

  //std::cout << "max_key: " << max_key << " | min_key: " << min_key << " | size: " << size << " | inverse_density: " << inverse_density << " | error_bound: " << temp_error << " | Q_value " << Q_value << std::endl;
  std::cout << " | size: " << size << " | inverse_density: " << inverse_density << " | build_time " << build_time << " | load_time " << weighted_sum <<" | loged_error " << adgMod::getQTableManagerInstance().Q_table[inverse_density].error_bound << " | error_bound: " << temp_error << " | meta's_err " << this->meta->error << " | Q_value " << Q_value << std::endl;
  adgMod::getQTableManagerInstance().Q_table[inverse_density].error_bound = temp_error;
  


  if (segs.empty()) return false;
  // fill in a dummy last segment (used in segment binary search)
  segs.push_back((Segment){temp, 0, 0});
  string_segments = std::move(segs);

  for (auto& str : string_segments) {
    // printf("%s %f\n", str.first.c_str(), str.second);
  }

  learned.store(true);
  // string_keys.clear();

  



  return true;
}

// static learning function to be used with LevelDB background scheduling
// level learning
void LearnedIndexData::LevelLearn(void* arg, bool nolock) {
  Stats* instance = Stats::GetInstance();
  bool success = false;
  bool entered = false;
  instance->StartTimer(8);

  VersionAndSelf* vas = reinterpret_cast<VersionAndSelf*>(arg);
  LearnedIndexData* self = vas->self;
  self->is_level = true;
  self->level = vas->level;
  Version* c;
  if (!nolock) {
    c = db->GetCurrentVersion();
  }
  if (db->version_count == vas->v_count) {
    entered = true;
    if (vas->version->FillLevel(adgMod::read_options, vas->level)) {
      self->filled = true;
      if (db->version_count == vas->v_count) {
        if (env->compaction_awaiting.load() == 0 && self->Learn()) {
          success = true;
        } else {
          self->learning.store(false);
        }
      }
    }
  }
  if (!nolock) {
    adgMod::db->ReturnCurrentVersion(c);
  }

  auto time = instance->PauseTimer(8, true);

  if (entered) {
    self->cost = time.second - time.first;
    learn_counter_mutex.Lock();
    events[1].push_back(new LearnEvent(time, 0, self->level, success));
    levelled_counters[6].Increment(vas->level, time.second - time.first);
    learn_counter_mutex.Unlock();
  }

  delete vas;
}

// static learning function to be used with LevelDB background scheduling
// file learning
uint64_t LearnedIndexData::FileLearn(void* arg) {
  Stats* instance = Stats::GetInstance();
  bool entered = false;
  instance->StartTimer(11);

  MetaAndSelf* mas = reinterpret_cast<MetaAndSelf*>(arg);
  LearnedIndexData* self = mas->self;
  self->level = mas->level;
  self->meta = mas->meta;

  Version* c = db->GetCurrentVersion();
  if (self->FillData(c, mas->meta)) {
    self->Learn();
    entered = true;
  } else {
    self->learning.store(false);
  }
  adgMod::db->ReturnCurrentVersion(c);

  auto time = instance->PauseTimer(11, true);

  if (entered) {
    // count how many file learning are done.
    self->cost = time.second - time.first;
    learn_counter_mutex.Lock();
    events[1].push_back(new LearnEvent(time, 1, self->level, true));
    levelled_counters[11].Increment(mas->level, time.second - time.first);
    learn_counter_mutex.Unlock();
  }

  if (!fresh_write) {
             self->WriteModel(adgMod::db->versions_->dbname_ + "/" +
             to_string(mas->meta->number) + ".fmodel");
             self->string_keys.clear();
             self->num_entries_accumulated.array.clear();
         }

  if (!fresh_write) delete mas->meta;
  delete mas;
  return entered ? time.second - time.first : 0;
}

// general model checker
bool LearnedIndexData::Learned() {
  if (learned_not_atomic)
    return true;
  else if (learned.load()) {
    learned_not_atomic = true;
    return true;
  } else
    return false;
}

// level model checker, used to be also learning trigger
bool LearnedIndexData::Learned(Version* version, int v_count, int level) {
  if (learned_not_atomic)
    return true;
  else if (learned.load()) {
    learned_not_atomic = true;
    return true;
  }
  return false;
  //        } else {
  //            if (level_learning_enabled && ++current_seek >= allowed_seek &&
  //            !learning.exchange(true)) {
  //                env->ScheduleLearning(&LearnedIndexData::Learn, new
  //                VersionAndSelf{version, v_count, this, level}, 0);
  //            }
  //            return false;
  //        }
}

// file model checker, used to be also learning trigger
bool LearnedIndexData::Learned(Version* version, int v_count,
                               FileMetaData* meta, int level) {
  if (learned_not_atomic)
    return true;
  else if (learned.load()) {
    learned_not_atomic = true;
    return true;
  } else
    return false;
  //        } else {
  //            if (file_learning_enabled && (true || level != 0 && level != 1)
  //            && ++current_seek >= allowed_seek && !learning.exchange(true)) {
  //                env->ScheduleLearning(&LearnedIndexData::FileLearn, new
  //                MetaAndSelf{version, v_count, meta, this, level}, 0);
  //            }
  //            return false;
  //        }
}

bool LearnedIndexData::FillData(Version* version, FileMetaData* meta) {
  // if (filled) return true;

  if (version->FillData(adgMod::read_options, meta, this)) {
    // filled = true;
    return true;
  }
  return false;
}

void LearnedIndexData::WriteModel(const string& filename) {
  if (!learned.load()) return;

  std::ofstream output_file(filename);
  output_file.precision(15);
  // output_file << adgMod::block_num_entries << " " << adgMod::block_size << " "
  //             << adgMod::entry_size << "\n";
  for (Segment& item : string_segments) {
    output_file << item.x << " " << item.k << " " << item.b << "\n";
  }
  output_file << 621;
  //             << " " << min_key << " " << max_key << " " << size << " " << level
  //             << " " << cost << "\n";
  // for (auto& pair : num_entries_accumulated.array) {
  //   output_file << pair.first << " " << pair.second << "\n";
  // }
}
void LearnedIndexData::ReadModel(const string& filename) {
    std::ifstream input_file(filename);

    if (!input_file.good()) return;

    std::string line;
    while (std::getline(input_file, line)) {
        std::istringstream iss(line);
        uint64_t x;
        double k = 0.0, b = 0.0;

        // 尝试读取 x, k, b
        if (!(iss >> x >> k >> b)) {
            // 检查是否遇到结束标记
            if (line == "621") {
                break;  // 读取到结束标记时停止读取
            } else {
                // std::cerr << "Error reading values from line: " << line << std::endl;
                continue;  // 跳过这一行，继续处理下一行
            }
        }

        // 将读取的数据添加到 string_segments_read 中
        string_segments.emplace_back(x, k, b);
    }

    learned.store(true);
}

void LearnedIndexData::ReportStats() {
  //        double neg_gain, pos_gain;
  //        if (num_neg_model == 0 || num_neg_baseline == 0) {
  //            neg_gain = 0;
  //        } else {
  //            neg_gain = ((double) time_neg_baseline / num_neg_baseline -
  //            (double) time_neg_model / num_neg_model) * num_neg_model;
  //        }
  //        if (num_pos_model == 0 || num_pos_baseline == 0) {
  //            pos_gain = 0;
  //        } else {
  //            pos_gain = ((double) time_pos_baseline / num_pos_baseline -
  //            (double) time_pos_model / num_pos_model) * num_pos_model;
  //        }

  printf("%d %d %lu %lu %lu\n", level, served, string_segments.size(), cost,
         size);  //, file_size);
  //        printf("\tPredicted: %lu %lu %lu %lu %d %d %d %d %d %lf\n",
  //        time_neg_baseline_p, time_neg_model_p, time_pos_baseline_p,
  //        time_pos_model_p,
  //                num_neg_baseline_p, num_neg_model_p, num_pos_baseline_p,
  //                num_pos_model_p, num_files_p, gain_p);
  //        printf("\tActual: %lu %lu %lu %lu %d %d %d %d %f\n",
  //        time_neg_baseline, time_neg_model, time_pos_baseline,
  //        time_pos_model,
  //               num_neg_baseline, num_neg_model, num_pos_baseline,
  //               num_pos_model, pos_gain + neg_gain);
}

void LearnedIndexData::FillCBAStat(bool positive, bool model, uint64_t time) {
  //        int& num_to_update = positive ? (model ? num_pos_model :
  //        num_pos_baseline) : (model ? num_neg_model : num_neg_baseline);
  //        uint64_t& time_to_update =  positive ? (model ? time_pos_model :
  //        time_pos_baseline) : (model ? time_neg_model : time_neg_baseline);
  //        time_to_update += time;
  //        num_to_update += 1;
}

LearnedIndexData* FileLearnedIndexData::GetModel(int number) {
  leveldb::MutexLock l(&mutex);
  if (file_learned_index_data.size() <= number)
    file_learned_index_data.resize(number + 1, nullptr);
  if (file_learned_index_data[number] == nullptr)
    file_learned_index_data[number] = new LearnedIndexData(file_allowed_seek, false);
  return file_learned_index_data[number];
}

bool FileLearnedIndexData::FillData(Version* version, FileMetaData* meta) {
  LearnedIndexData* model = GetModel(meta->number);
  return model->FillData(version, meta);
}

std::vector<std::string>& FileLearnedIndexData::GetData(FileMetaData* meta) {
  auto* model = GetModel(meta->number);
  return model->string_keys;
}

bool FileLearnedIndexData::Learned(Version* version, FileMetaData* meta,
                                   int level) {
  LearnedIndexData* model = GetModel(meta->number);
  return model->Learned(version, db->version_count, meta, level);
}

AccumulatedNumEntriesArray* FileLearnedIndexData::GetAccumulatedArray(
    int file_num) {
  auto* model = GetModel(file_num);
  return &model->num_entries_accumulated;
}

std::pair<uint64_t, uint64_t> FileLearnedIndexData::GetPosition(
    const Slice& key, int file_num) {
  return file_learned_index_data[file_num]->GetPosition(key);
}

FileLearnedIndexData::~FileLearnedIndexData() {
  leveldb::MutexLock l(&mutex);
  for (auto pointer : file_learned_index_data) {
    delete pointer;
  }
}

void FileLearnedIndexData::Report() {
  leveldb::MutexLock l(&mutex);

  std::set<uint64_t> live_files;
  adgMod::db->versions_->AddLiveFiles(&live_files);

  for (size_t i = 0; i < file_learned_index_data.size(); ++i) {
    auto pointer = file_learned_index_data[i];
    if (pointer != nullptr && pointer->cost != 0) {
      printf("FileModel %lu %d ", i, i > watermark);
      pointer->ReportStats();
    }
  }
}

void AccumulatedNumEntriesArray::Add(uint64_t num_entries, string&& key) {
  array.emplace_back(num_entries, key);
}

bool AccumulatedNumEntriesArray::Search(const Slice& key, uint64_t lower,
                                        uint64_t upper, size_t* index,
                                        uint64_t* relative_lower,
                                        uint64_t* relative_upper) {
  if (adgMod::MOD == 4) {
    uint64_t lower_pos = lower / array[0].first;
    uint64_t upper_pos = upper / array[0].first;
    if (lower_pos != upper_pos) {
      while (true) {
        if (lower_pos >= array.size()) return false;
        if (key <= array[lower_pos].second) break;
        lower = array[lower_pos].first;
        ++lower_pos;
      }
      upper = std::min(upper, array[lower_pos].first - 1);
      *index = lower_pos;
      *relative_lower =
          lower_pos > 0 ? lower - array[lower_pos - 1].first : lower;
      *relative_upper =
          lower_pos > 0 ? upper - array[lower_pos - 1].first : upper;
      return true;
    }
    *index = lower_pos;
    *relative_lower = lower % array[0].first;
    *relative_upper = upper % array[0].first;
    return true;

  } else {
    size_t left = 0, right = array.size() - 1;
    while (left < right) {
      size_t mid = (left + right) / 2;
      if (lower < array[mid].first)
        right = mid;
      else
        left = mid + 1;
    }

    if (upper >= array[left].first) {
      while (true) {
        if (left >= array.size()) return false;
        if (key <= array[left].second) break;
        lower = array[left].first;
        ++left;
      }
      upper = std::min(upper, array[left].first - 1);
    }

    *index = left;
    *relative_lower = left > 0 ? lower - array[left - 1].first : lower;
    *relative_upper = left > 0 ? upper - array[left - 1].first : upper;
    return true;
  }
}

bool AccumulatedNumEntriesArray::SearchNoError(uint64_t position, size_t* index,
                                               uint64_t* relative_position) {
  *index = position / array[0].first;
  *relative_position = position % array[0].first;
  return *index < array.size();

  //        size_t left = 0, right = array.size() - 1;
  //        while (left < right) {
  //            size_t mid = (left + right) / 2;
  //            if (position < array[mid].first) right = mid;
  //            else left = mid + 1;
  //        }
  //        *index = left;
  //        *relative_position = left > 0 ? position - array[left - 1].first :
  //        position; return left < array.size();
}

uint64_t AccumulatedNumEntriesArray::NumEntries() const {
  return array.empty() ? 0 : array.back().first;
}

}  // namespace adgMod
