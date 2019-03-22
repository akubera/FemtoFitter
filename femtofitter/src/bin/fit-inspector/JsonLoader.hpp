///
/// \file fit-inspector/JsonLoader.hpp
///

#pragma once

#ifndef FITINSPECTOR_JSONLOADER_HPP
#define FITINSPECTOR_JSONLOADER_HPP

#include <TString.h>

// #include <variant>


struct JsonLoadOk {
  TString data;
};

struct JsonLoadError {
  TString reason;
};

// using JsonLoadResult = std::variant<JsonLoadOk, JsonLoadError>;

struct JsonLoadResult {
  bool is_error;
  union {
    JsonLoadOk ok;
    JsonLoadError err;
  };

  static JsonLoadResult Error(TString msg)
    {
      return JsonLoadResult(true, msg);
    }

  static JsonLoadResult Ok(TString msg)
    {
      return JsonLoadResult(false, msg);
    }

  JsonLoadResult(const JsonLoadResult &o)
    : is_error(o.is_error)
    {
      if (is_error) {
        err.reason = o.err.reason;
      } else {
        ok.data = o.ok.data;
      }
    }

  JsonLoadResult& operator=(const JsonLoadResult &o) = delete;

  ~JsonLoadResult()
    {
      // if (is_error) {
      //   delete &err;
      // } else {
      //   delete &ok;
      // }
    }

  operator bool() const
    { return !is_error; }

protected:
  JsonLoadResult(bool is_err, TString msg)
    : is_error(is_err)
    , ok({msg})
    {
      if (is_error) {
        delete &ok;
      //   // err.reason = msg;
        new (&err) JsonLoadError { msg };
      }
      // } else {
      //   new (&ok) JsonLoadOk {msg};
      // }
    }
};

#endif
