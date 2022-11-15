#!/bin/bash

proc="$1"
particle="$2"

root -l -b -q "AnaClusters.C(${proc}, \"${particle}\")"
