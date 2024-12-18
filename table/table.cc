// Copyright (c) 2011 The LevelDB Authors. All rights reserved.
// Use of this source code is governed by a BSD-style license that can be
// found in the LICENSE file. See the AUTHORS file for names of contributors.

#include <iostream>
#include "leveldb/table.h"

#include "leveldb/cache.h"
#include "leveldb/comparator.h"
#include "leveldb/env.h"
#include "leveldb/filter_policy.h"
#include "leveldb/options.h"
#include "table/block.h"
#include "table/filter_block.h"
#include "table/format.h"
#include "table/two_level_iterator.h"
#include "util/coding.h"
#include "mod/stats.h"
#include "mod/learned_index.h"
#include <fstream>
#include "../db/version_set.h"

namespace leveldb {

Table::Rep::~Rep() {
    delete filter;
    delete[] filter_data;
    delete index_block;
}

Status Table::Open(const Options& options, RandomAccessFile* file,
                   uint64_t size, Table** table) {
  *table = nullptr;
  if (size < Footer::kEncodedLength) {
    return Status::Corruption("file is too short to be an sstable");
  }

  char footer_space[Footer::kEncodedLength];
  Slice footer_input;
  Status s = file->Read(size - Footer::kEncodedLength, Footer::kEncodedLength,
                        &footer_input, footer_space);
  if (!s.ok()) return s;

  Footer footer;
  s = footer.DecodeFrom(&footer_input);
  if (!s.ok()) return s;

  // Read the index block
  BlockContents index_block_contents;
  if (s.ok()) {
    ReadOptions opt;
    if (options.paranoid_checks) {
      opt.verify_checksums = true;
    }
    s = ReadBlock(file, opt, footer.index_handle(), &index_block_contents);
  }

  if (s.ok()) {
    // We've successfully read the footer and the index block: we're
    // ready to serve requests.
    Block* index_block = new Block(index_block_contents);
    Rep* rep = new Table::Rep;
    rep->options = options;
    rep->file = file;
    rep->metaindex_handle = footer.metaindex_handle();
    rep->index_block = index_block;
    rep->cache_id = (options.block_cache ? options.block_cache->NewId() : 0);
    rep->filter_data = nullptr;
    rep->filter = nullptr;
    *table = new Table(rep);
    (*table)->ReadMeta(footer);
  }

  return s;
}

void Table::ReadMeta(const Footer& footer) {
  if (rep_->options.filter_policy == nullptr) {
    return;  // Do not need any metadata
  }

  // TODO(sanjay): Skip this if footer.metaindex_handle() size indicates
  // it is an empty block.
  ReadOptions opt;
  if (rep_->options.paranoid_checks) {
    opt.verify_checksums = true;
  }
  BlockContents contents;
  if (!ReadBlock(rep_->file, opt, footer.metaindex_handle(), &contents).ok()) {
    // Do not propagate errors since meta info is not needed for operation
    return;
  }
  Block* meta = new Block(contents);

  Iterator* iter = meta->NewIterator(BytewiseComparator());
  std::string key = "filter.";
  key.append(rep_->options.filter_policy->Name());
  iter->Seek(key);
  if (iter->Valid() && iter->key() == Slice(key)) {
    ReadFilter(iter->value());
  }
  delete iter;
  delete meta;
}

void Table::ReadFilter(const Slice& filter_handle_value) {
  Slice v = filter_handle_value;
  BlockHandle filter_handle;
  if (!filter_handle.DecodeFrom(&v).ok()) {
    return;
  }

  // We might want to unify with ReadBlock() if we start
  // requiring checksum verification in Table::Open.
  ReadOptions opt;
  if (rep_->options.paranoid_checks) {
    opt.verify_checksums = true;
  }
  BlockContents block;
  if (!ReadBlock(rep_->file, opt, filter_handle, &block).ok()) {
    return;
  }
  if (block.heap_allocated) {
    rep_->filter_data = block.data.data();  // Will need to delete later
  }
  rep_->filter = new FilterBlockReader(rep_->options.filter_policy, block.data);
}

Table::~Table() { delete rep_; }

static void DeleteBlock(void* arg, void* ignored) {
  delete reinterpret_cast<Block*>(arg);
}

static void DeleteCachedBlock(const Slice& key, void* value) {
  Block* block = reinterpret_cast<Block*>(value);
  delete block;
}

static void ReleaseBlock(void* arg, void* h) {
  Cache* cache = reinterpret_cast<Cache*>(arg);
  Cache::Handle* handle = reinterpret_cast<Cache::Handle*>(h);
  cache->Release(handle);
}

// Convert an index iterator value (i.e., an encoded BlockHandle)
// into an iterator over the contents of the corresponding block.
Iterator* Table::BlockReader(void* arg, const ReadOptions& options,
                             const Slice& index_value) {
  Table* table = reinterpret_cast<Table*>(arg);
  Cache* block_cache = table->rep_->options.block_cache;
  Block* block = nullptr;
  Cache::Handle* cache_handle = nullptr;

  BlockHandle handle;
  Slice input = index_value;
  Status s = handle.DecodeFrom(&input);
  // We intentionally allow extra stuff in index_value so that we
  // can add more features in the future.

  if (s.ok()) {
    BlockContents contents;
    // if (block_cache != nullptr) {
    //   char cache_key_buffer[16];
    //   EncodeFixed64(cache_key_buffer, table->rep_->cache_id);
    //   EncodeFixed64(cache_key_buffer + 8, handle.offset());
    //   Slice key(cache_key_buffer, sizeof(cache_key_buffer));
    //   cache_handle = block_cache->Lookup(key);
    //   if (cache_handle != nullptr) {
    //     block = reinterpret_cast<Block*>(block_cache->Value(cache_handle));
    //   } else {
    //     s = ReadBlock(table->rep_->file, options, handle, &contents);
    //     if (s.ok()) {
    //       block = new Block(contents);
    //       if (contents.cachable && options.fill_cache) {
    //         cache_handle = block_cache->Insert(key, block, block->size(),
    //                                            &DeleteCachedBlock);
    //       }
    //     }
    //   }
    // } else {
      s = ReadBlock(table->rep_->file, options, handle, &contents);
      if (s.ok()) {
        block = new Block(contents);
      }
    // }
  }

  Iterator* iter;
  if (block != nullptr) {
    iter = block->NewIterator(table->rep_->options.comparator);
    if (cache_handle == nullptr) {
      iter->RegisterCleanup(&DeleteBlock, block, nullptr);
    } else {
      iter->RegisterCleanup(&ReleaseBlock, block_cache, cache_handle);
    }
  } else {
    iter = NewErrorIterator(s);
  }
  return iter;
}

Iterator* Table::NewIterator(const ReadOptions& options) const {
  return NewTwoLevelIterator(
      rep_->index_block->NewIterator(rep_->options.comparator),
      &Table::BlockReader, const_cast<Table*>(this), options);
}
Status Table::InternalGet(const ReadOptions& options, const Slice& k, void* arg,
                          void (*handle_result)(void*, const Slice&, const Slice&), int level,
                          FileMetaData* meta, uint64_t lower, uint64_t upper, bool learned, Version* version) {
  adgMod::Stats* instance = adgMod::Stats::GetInstance();
  Status s;
  Iterator* iiter = rep_->index_block->NewIterator(rep_->options.comparator);
  ParsedInternalKey parsed_key;
  ParseInternalKey(k, &parsed_key);


#ifdef INTERNAL_TIMER
    instance->StartTimer(2);
#endif
    iiter->Seek(k);
#ifdef INTERNAL_TIMER
    instance->PauseTimer(2);
#endif

    if (iiter->Valid()) {
#ifdef INTERNAL_TIMER
      instance->StartTimer(15);
#endif
      Slice handle_value = iiter->value();
      FilterBlockReader* filter = rep_->filter;
      BlockHandle handle;
      if (filter != nullptr && handle.DecodeFrom(&handle_value).ok() &&
          !filter->KeyMayMatch(handle.offset(), k)) {

#ifdef INTERNAL_TIMER
        auto time = instance->PauseTimer(15, true);
        adgMod::levelled_counters[9].Increment(level, time.second - time.first);
#endif
        // Not found
      } else {
#ifdef INTERNAL_TIMER
        auto time = instance->PauseTimer(15, true);
        adgMod::levelled_counters[9].Increment(level, time.second - time.first);
        instance->StartTimer(5);
#endif
        Iterator* block_iter = BlockReader(this, options, iiter->value());
#ifdef INTERNAL_TIMER
        instance->PauseTimer(5);
        instance->StartTimer(3);
#endif
        block_iter->Seek(k);
#ifdef INTERNAL_TIMER
        instance->PauseTimer(3);
#endif
        if (block_iter->Valid()) {
          (*handle_result)(arg, block_iter->key(), block_iter->value());
        }
        s = block_iter->status();
        delete block_iter;
      }
    }
    if (s.ok()) {
      
      s = iiter->status();
    }
    delete iiter;
    return s;
}

uint64_t Table::ApproximateOffsetOf(const Slice& key) const {
  Iterator* index_iter =
      rep_->index_block->NewIterator(rep_->options.comparator);
  index_iter->Seek(key);
  uint64_t result;
  if (index_iter->Valid()) {
    BlockHandle handle;
    Slice input = index_iter->value();
    Status s = handle.DecodeFrom(&input);
    if (s.ok()) {
      result = handle.offset();
    } else {
      // Strange: we can't decode the block handle in the index block.
      // We'll just return the offset of the metaindex block, which is
      // close to the whole file size for this case.
      result = rep_->metaindex_handle.offset();
    }
  } else {
    // key is past the last key in the file.  Approximate the offset
    // by returning the offset of the metaindex block (which is
    // right near the end of the file).
    result = rep_->metaindex_handle.offset();
  }
  delete index_iter;
  return result;
}





void Table::FillData(const ReadOptions& options, adgMod::LearnedIndexData* data) {
    if (data->filled) return;

    // 출력 파일 스트림 객체 생성

    
    // std::ofstream output_file("output_keys.txt", std::ios::out | std::ios::app);
    // if (!output_file.is_open()) {
    //     // 파일 열기 실패 시 적절히 처리
    //     throw std::runtime_error("Failed to open output_keys.txt");
    // }

    Status status;
    Block::Iter* index_iter = dynamic_cast<Block::Iter*>(rep_->index_block->NewIterator(rep_->options.comparator));

    for (uint32_t i = 0; i < index_iter->num_restarts_; ++i) {
        index_iter->SeekToRestartPoint(i);
        index_iter->ParseNextKey();
        assert(index_iter->Valid());
        Block::Iter* block_iter = dynamic_cast<Block::Iter*>(BlockReader(this, options, index_iter->value()));

        ParsedInternalKey parsed_key;
        int num_entries_this_block = 0;
        for (block_iter->SeekToRestartPoint(0); block_iter->ParseNextKey(); ++num_entries_this_block) {
            ParseInternalKey(block_iter->key(), &parsed_key);
            data->string_keys.emplace_back(parsed_key.user_key.data(), parsed_key.user_key.size());

            // 파일에 user_key 기록
            // output_file << std::string(parsed_key.user_key.data()) << "\n";
        }

        if (!adgMod::block_num_entries_recorded) {
            adgMod::block_num_entries = num_entries_this_block;
            adgMod::block_num_entries_recorded = true;
            adgMod::entry_size = block_iter->restarts_ / num_entries_this_block;
            BlockHandle temp;
            Slice temp_slice = index_iter->value();
            temp.DecodeFrom(&temp_slice);
            adgMod::block_size = temp.size() + kBlockTrailerSize;
        }

        uint64_t current_total = data->num_entries_accumulated.NumEntries();
        if (!data->Learned()) {
            data->num_entries_accumulated.Add(current_total + num_entries_this_block, std::string(parsed_key.user_key.data(), parsed_key.user_key.size()));
        }
        delete block_iter;
    }

    data->filled = true;
    delete index_iter;

    // 파일 스트림 닫기
    // output_file.close();
}


}  // namespace leveldb
