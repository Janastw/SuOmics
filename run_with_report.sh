#!/bin/bash

mkdir -p results/logs
nextflow run main.nf \
  -with-report results/logs/report.html \
  -with-trace results/logs/trace.txt \
  -with-timeline results/logs/timeline.html \
  -with-dag results/logs/flowchart.html