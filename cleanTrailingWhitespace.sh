#!/bin/bash

find Derivs -type f -print0 | xargs -0 sed -i .bak -E "s/[[:space:]]*$//"
find Matlab -type f -print0 | xargs -0 sed -i .bak -E "s/[[:space:]]*$//"
