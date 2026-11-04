#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Apr 12 16:36:53 2020
"""

import bundleTools as bt
import bundleTools3 as bt3
import os
import sys
import numpy as np

if len(sys.argv) != 4:
 	sys.exit('Argument error. Correct syntax: "python ' + sys.argv[0] + ' input_bundle.bundles output_bundle.bundles min_size_in_mm"');

raw_bundle_path = sys.argv[1];
filtered_bundle_path = sys.argv[2];
min_size = int(sys.argv[3]);

raw_bundle = bt.read_bundle(raw_bundle_path);

filtered_bundle = [];

for fiber in raw_bundle:
    a = fiber[0];
    b = fiber[1];
    dist = np.linalg.norm(a-b);
    length = dist*20;
    
    if length > min_size:
        filtered_bundle.append(fiber);
        
bt3.write_bundle(filtered_bundle_path, filtered_bundle);
