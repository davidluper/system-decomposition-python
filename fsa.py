#############################################################################################
#System Decomposition
#by davidluper
#
#Copyright 2010 davidluper
#This file is part of System Decomposition.
#
#System Decomposition is free software: you can redistribute it and/or modify it under the
#terms of the GNU General Public License as published by the Free Software Foundation,
#either version 3 of the License, or (at your option) any later version.
#
#System Decomposition is distributed in the hope that it will be useful, but WITHOUT ANY
#WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
#PURPOSE. See the GNU General Public License for more details.
#
#You should have received a copy of the GNU General Public License along with
#System Decomposition.  If not, see http://www.gnu.org/licenses/.
#############################################################################################

import sys;
import os;
import shutil;
from numpy import *;
import system_decomp;

def Normalize(l,multiplier):
    total = 0.0;
    for i in l: total += float(i);
    for i in range(len(l)): l[i] = float(l[i]) / float(total);
    for i in range(len(l)): l[i] = float(multiplier) * float(l[i]);
    return l;

def PathCount(_dir):
    cnt = int(0);

    _file = None;
    print(_dir);
    for item in os.listdir(_dir):
        if "path" in item and ".dat" in item:
            _file = item;
            break;

    if _file is None:
        return None;

    _file = _dir + _file;
    
    pathfile = open(_file, 'r');
    for line in pathfile:
        cnt += 1;
    pathfile.close();

    print(_file + "\t" + str(cnt));
    return float(cnt);

def MakeDistanceFile(dir1, dir2, multiplier1, multiplier2, inline):
    hist1 = matrix(Normalize(system_decomp.LoadHistogram(dir1 + "freqcount.txt"), multiplier1));
    hist2 = matrix(Normalize(system_decomp.LoadHistogram(dir2 + "freqcount.txt"), multiplier2));

    if inline:
        dvm = system_decomp.LoadInlineDVM(dir1 + "dvm-inlinesearch.txt");
    else:
        dvm = system_decomp.LoadDVM(dir1 + "dv-log.txt", dir1 + "freqcount_flux.txt");

    #hist2mod is a modified version of hist2 where hist2mod is the corresponding point to hist1 with the minimum distance
    dist, hist2mod = system_decomp.DVMDistanceCalc(hist1, hist2, dvm);

    fluxes = [];
    fluxfile = open(dir1 + "freqcount_flux.txt", 'r');
    for line in fluxfile:
        fluxes.append(line.strip().replace(",", "_"));
    fluxfile.close();

    label = "absolute";
    if float(multiplier1) == float(1.0) and float(multiplier2) == float(1.0):
        label = "relative";
    outfile = open(dir1 + "distances_" + label + ".csv", 'w');
    outfile.write("System 1 Point,System 2 Point,Flux,Distance\n");
    for i in range(len(fluxes)):
        outfile.write(str(hist1[0,i]) + "," + str(hist2mod[0,i]) + "," + fluxes[i] + "," + str(math.fabs(float(hist1[0,i]) - float(hist2mod[0,i]))) + "\n");

    outfile.close();
    return;

#sysargs
# 1 - [dir1 (must have freqcount.txt, freqcount_flux.txt, and dv-log.txt)]
# 2 - [dir2 (must have freqcount.txt from second run)]

folderdelim = ("/", "\\")[sys.argv[1].find("\\") > -1];

dir1 = (sys.argv[1] + folderdelim, sys.argv[1])[sys.argv[1].endswith(folderdelim)];
dir2 = (sys.argv[2] + folderdelim, sys.argv[2])[sys.argv[2].endswith(folderdelim)];

cnt1 = PathCount(dir1);
cnt2 = PathCount(dir2);

if cnt1 is None:
    print("dir1 does not contain a paths.dat file");
    sys.exit(0);
    
if cnt2 is None:
    print("dir2 does not contain a paths.dat file");
    sys.exit(0);

useinline = bool(0);

if len(sys.argv) == 4 and sys.argv[3] == "-inline":
    useinline = bool(1);
    
MakeDistanceFile(dir1, dir2, float(1.0), float(1.0), useinline);
MakeDistanceFile(dir1, dir2, cnt1, cnt2, useinline);
