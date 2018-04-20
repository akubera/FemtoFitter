///
/// \file femtofitter.h
///

#pragma once

#ifndef FEMTOFITTER_H
#define FEMTOFITTER_H


// ====== C-Macros

#define LOOP_BINS_3D_HIST(START, STOP, h)           \
  for (uint_fast16_t k=START; k<=STOP; ++k) {       \
    double ql = h.GetZaxis()->GetBinCenter(k);      \
    for (uint_fast16_t j=START; j<=STOP; ++j) {     \
      double qs = h.GetYaxis()->GetBinCenter(j);    \
      for (uint_fast16_t i=START; i<=STOP; ++i) {   \
        double qo = h.GetXaxis()->GetBinCenter(i);


#define LOOP_BINS_3D_VEC(START, STOP, obj)        \
  for (uint_fast16_t k=START; k<=STOP; ++k) {     \
    double ql = q_vec[k];                         \
    for (uint_fast16_t j=START; j<=STOP; ++j) {   \
      double qs = q_vec[j];                       \
      for (uint_fast16_t i=START; i<=STOP; ++i) { \
        double qo = q_vec[i];


#define END_LOOP_BINS_3D } } }


// ======== include

#include "femtomath.h"
#include """

#endif
