///
/// \file FileManager.hpp
///

#pragma once

#include <map>
#include <string>
#include <memory>

#include <TFile.h>


struct FileManager {

  std::map<std::string, std::shared_ptr<TFile>> filemap;

  std::shared_ptr<TFile> current_file;


  void add(std::string fname, std::shared_ptr<TFile> ptr)
    {
       filemap[fname] = ptr;
       current_file = ptr;
    }

  std::shared_ptr<TFile> get(std::string fname)
    {
      auto found = filemap.find(fname);
      if (found == filemap.end()) {
        return nullptr;
      }
      return found->second;
    }


};
