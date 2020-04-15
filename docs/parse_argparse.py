#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
COPYRIGHT
Copyright 2020 Ian J. Miller, Evan R. Rees, Kyle Wolf, Siddharth Uppal,
Shaurya Chanana, Izaak Miller, Jason C. Kwan

This file is part of Autometa.

Autometa is free software: you can redistribute it and/or modify
it under the terms of the GNU Affero General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

Autometa is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU Affero General Public License for more details.

You should have received a copy of the GNU Affero General Public License
along with Autometa. If not, see <http://www.gnu.org/licenses/>.
COPYRIGHT

Script used to parse the arparse block of all the autometa scripts, copy them to a new
file and run the help text (--help) from there. 
Used as an alternative to sphinxcontrib-programoutput which used all the heavy dependencies
preventing the intergration with readthedocs.
"""


import glob
import subprocess
import os
import textwrap 
import shutil


path_scripts = '../docs/source/scripts'
if os.path.exists(path_scripts):
    shutil.rmtree(path_scripts) 

for dir_name in glob.glob("../autometa/**/*.py", recursive=True): 
    if "__" in dir_name: #to not include __init__ files
        continue
    print(dir_name)
    with open (dir_name) as fh:
        path_temp_write = '../docs/temp_write.py'
        fh_write = open(path_temp_write, 'w+')            
        for line in fh:
            if line.strip() == 'import argparse':
                line = line.strip()
                fh_write.write(line + '\n')
                break
        # Reads text until the end of the block:
        for line in fh:  # This keeps reading the file
            if line.strip() == 'args = parser.parse_args()':
                line = line.strip()
                fh_write.write(line + '\n')
                break
            line = line.strip()
            fh_write.write(line + '\n')
        fh_write.close()
    
    basename = os.path.basename(dir_name) #Extract the filename, eg: kmers.py
    #print (basename)
    len_base = len(basename) #needed to get the number of '=' in the header
    file_name_base = os.path.splitext(basename)[0] #extract the kmers of kmers.py
    file_name = file_name_base + ".rst" # make it kmers.rst   
    #print(file_name)
    main_path = os.path.dirname(os.path.abspath(dir_name)) # extract '/home/siddharth/Autometa/autometa/common/external'
    #print(main_path)
    path = main_path.split('/')
    #print(main_path)
    main_path = path [-1] # extraxt 'external'
    path[path.index('autometa')] = 'docs/source/scripts' #replace autometa with the specified path in the list
    path_ori=path='/'.join(path)

    if not os.path.exists(path): #path now has 'docs/source/scripts/comon/external' instead of autometa
        #make a new directory as per the structure in autometa
        os.makedirs(path) #take care now, this will make each folder start with a lower case letter
    path = os.path.join(path, file_name) #adding the kmers.rst to the above complete path
    #now the path is 'docs/source/scripts/comon/external/samtools.rst'
    #print(path)

    cmd = 'python3 '+  os.path.abspath(path_temp_write) + ' -h'
    #print(cmd)

    path_temp_write_rst = os.path.abspath(path_temp_write)
    #print(path_temp_write_rst)
    base_path_temp_write = os.path.splitext(os.path.abspath(path_temp_write))[0]
    #print(base_path_temp_write)
    path_temp_write_rst = base_path_temp_write + '.rst'
    #print(path_temp_write_rst)

    #this block writes the argparse output to the respective '.rst' file with proper indentation
    #this is the output that will be under '..code-block::' and will be displayed in html
    with open(path_temp_write_rst, 'w+') as stdout:
        subprocess.call(cmd, stdout=stdout, shell=True) #capture the argparse output
        stdout.seek(0) #takes the cursor to beginning. Copy and indent stdout.
        #The indentation is done because the arparse text needs to be indented
        #under the ..code-block:: section in the final rst file
        wrapper = textwrap.TextWrapper(initial_indent='\t', subsequent_indent='\t', width = 700)
        wrapped=""
        with open(path,'w') as fh_rst:
            fh_rst.write('='*len_base + '\n' + basename + '\n' + '='*len_base + '\n' + '\n'+ '.. code-block:: shell' + '\n \n')
            for line in stdout:
                wrapped = wrapped + wrapper.fill(line) + '\n' #copying the indented text
            fh_rst.write(wrapped) #writing the indented tex to the final rst files, this will be dislayed in html
    
    #print(main_path)
    main_file_no_rst = main_path + '_main' #making a file external_main
    main_file = main_file_no_rst + '.rst' #making file external_main.rst
    #print(main_file)
    main_file = os.path.join(path_ori, main_file) # we now have /home/siddharth/Autometa/docs/source/scripts/common/external/external_main.rst
    #print(main_file)
    #print(file_name_base)

    #this block will make the file 'Usage.rst' this is the script that will be called by the toctree in index.rst
    if not os.path.exists('../docs/source/scripts/Usage.rst'):
        with open ('../docs/source/scripts/Usage.rst', 'w') as fh_main:
            fh_main.write('='*len('Usage') + '\n' + 'Usage \n' + '='*len('Usage') + '\n' + '\n'+ '.. toctree::' +
            '\n \t :maxdepth: 3' + '\n \t :caption: Table of Contents \n')

    #this will create a _main file for each directory. These rst files will be called by the toctreee in scripts_main.rst
    if not os.path.exists(main_file):
        len_head = len(main_path)
        with open (main_file, 'w+') as fh_main:
            fh_main.write('='*len_head + '\n' + main_path + '\n' + '='*len_head + '\n' + '\n'+ '.. toctree::' +
            '\n \t :maxdepth: 2' + '\n \t :caption: Table of Contents \n' )

        #need to create a link to external directory which is within the common directory
        if main_path == 'common':
            with open (main_file, 'a+') as fh_main:
                fh_main.write(f' \n \t external/external_main')

        #create links to all main scripts in the scripts_main toctress, eg: external/external_main.rst (in toctree)
        with open ('../docs/source/scripts/Usage.rst', 'a') as fh_main:
            if main_path == 'external':
                #Common will link the external directory we do not need to include it in the main scripts_main.rst
                pass
            else:
                fh_main.write( f'\n \t {main_path}/{main_file_no_rst}')

    #open each _main.rst and add the name of the scripts within it
    #eg in external.rst, we add samtool, bedtools, prodigal, etc
    with open (main_file, 'a+') as fh_main:
            fh_main.write(f'\n \t {file_name_base}')
        
os.remove(path_temp_write)
os.remove(base_path_temp_write+'.rst')