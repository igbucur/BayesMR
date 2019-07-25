#!/usr/bin/env bash
# Script for running experiments in SMMR paper

# 6.1 Birth weight versus fasting glucose
bin/polychord_MR ini/BW_FG_dir.ini ini/BW_FG_model_dir.ini
bin/polychord_MR ini/BW_FG_rev.ini ini/BW_FG_model_rev.ini

# 6.2 risk of Parkinson's disease


# 6.3 Caffeine vs Nicotine

