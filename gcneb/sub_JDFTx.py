#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov  6 20:43:05 2020

@author: zaba1157
"""
import argparse
import os
import subprocess

opj = os.path.join

try:
    modules=' '.join(os.environ['JDFTx_mods'].split('_'))
except:
    modules='' #'comp-intel/2020.1.217 intel-mpi/2020.1.217 cuda/10.2.89 vasp/6.1.1 mkl/2020.1.217 gsl/2.5/gcc openmpi/4.0.4/gcc-8.4.0 gcc/7.4.0'

try:
    comp=os.environ['JDFTx_Computer']
except:
    comp='Eagle'

def write(nodes,cores,time,out,alloc,qos,script,short_recursive,procs,gpu,testing):
#    if short_recursive == 'True':
#        if time != 4: 
#            print('Time limit set to 04:00:00 from short_recursive')
#            time = 4
#        if nodes > 4: 
#            print('Nodes set to 4 from short_recursive')
#            nodes = 4
#        out = 'sc_'+out
    testing = True if testing == 'True' else False
    np=nodes*cores
    writelines = '#!/bin/bash'+'\n'
    writelines+='#SBATCH -J '+out+'\n'
    if testing:
        if comp == 'Summit':
            writelines+='#SBATCH --time=0:30:00'+'\n'
        elif comp == 'Eagle':
            writelines+='#SBATCH --time=1:00:00'+'\n'
    else:
        writelines+='#SBATCH --time='+str(time)+':00:00'+'\n'
    writelines+='#SBATCH -o '+out+'-%j.out'+'\n'
    writelines+='#SBATCH -e '+out+'-%j.err'+'\n'
    
    if comp == 'Eagle' and gpu != 'True':
        writelines+='#SBATCH -N 1 -n 2 -c 18 --hint=nomultithread'+'\n'
    else:
        writelines+='#SBATCH --tasks '+str(np)+'\n'
        writelines+='#SBATCH --nodes '+str(nodes)+'\n'
        writelines+='#SBATCH --ntasks-per-node '+str(cores)+'\n'

    if alloc=='environ':
        writelines+='#SBATCH --account='+os.environ['JDFTx_allocation']+'\n'
    else:
        writelines+='#SBATCH --account='+alloc+'\n'
    if qos=='high' and comp == 'Eagle':
        writelines+='#SBATCH --qos=high'+'\n'

    if (time == 1 or testing) and comp == 'Eagle':
        writelines+='#SBATCH --partition=debug\n'
        
    if comp == 'Summit':
        if gpu == 'True':
            writelines+='#SBATCH --partition sgpu\n'
        elif testing:
            writelines+='#SBATCH --partition shas-testing\n'
        else:
            writelines+='#SBATCH --partition shas\n'    
    
#    if comp == 'Eagle' and gpu != 'True':
#        writelines+='#SBATCH --hint=nomultithread'
    
    if comp != 'Eagle':
        writelines+='\nexport JDFTx_NUM_PROCS='+str(procs)+'\n' # previously np
    if comp == 'Summit':
        writelines+='SLURM_EXPORT_ENV=ALL\n'

    #writelines+='module load comp-intel/2020.1.217 intel-mpi/2020.1.217 cuda/10.2.89 vasp/6.1.1 mkl/2020.1.217 gsl/2.5/gcc openmpi/4.0.4/gcc-8.4.0 gcc/7.4.0'+'\n\n'
    if comp == 'Eagle':
        writelines+='\n'+'module use -a /nopt/nrel/apps/modules/test/modulefiles'+'\n'
#        writelines+='module load gcc/8.4.0 openmpi/4.1.1/gcc+cuda hdf5/1.10.7/gcc-ompi gsl/2.5/gcc cmake mkl'
    
    if modules != '':
        writelines+='\nmodule load '+modules+'\n\n'

    if short_recursive == 'True': # removed time constraint
        # short_recursive command runs timer script before and after JDFT to check if walltime is hit
        writelines+='timeout 10 python '+os.path.join(jdftx_python_dir,'timer.py')+' > timer'+'\n'
        
        time_lim = str(3600 * time - 100)
        if gpu == 'True':
            writelines+='timeout '+time_lim+' ' + 'python '+script+' -g True > '+out+'\n'
        else:
            writelines+='timeout '+time_lim+' ' + 'python '+script+' > '+out+'\n'
        writelines+='timeout 10 python '+os.path.join(jdftx_python_dir,'timer.py')+' > timer'+'\n'

    else:
        if gpu == 'True':
            writelines+='python '+script+' -g True > '+out+'\n'
        else:
            writelines+='python '+script+' > '+out+'\n'
    writelines+='exit 0'+'\n'

    with open('submit.sh','w') as f:
        f.write(writelines)

# bridges requires --nodes, -t, --ntasks-per-node, -p
def write_bridges(nodes,cores,time,out,partition,qos,script,short_recursive,procs):
#    if short_recursive == 'True':
#        if time > 4: 
#            print('Time limit set to 04:00:00 from short_recursive')
#            time = 4
#        out = 'sc_'+out

    if partition == 'RM':
        cores = 128
    elif partition == 'RM-shared':
        cores = 64
#    np=nodes*cores
    
    writelines = '#!/bin/bash'+'\n'
    writelines+='#SBATCH -J '+out+'\n'
    writelines+='#SBATCH -t '+str(time)+':00:00'+'\n'
    writelines+='#SBATCH -o '+out+'-%j.out'+'\n'
    writelines+='#SBATCH -e '+out+'-%j.err'+'\n'
    writelines+='#SBATCH -p '+partition+'\n'
    writelines+='#SBATCH --nodes '+str(nodes)+'\n'
    writelines+='#SBATCH --ntasks-per-node '+str(cores)+'\n'

    if qos=='high':
        writelines+='#SBATCH --qos=high'+'\n'
    
    writelines+='\nexport JDFTx_NUM_PROCS='+str(procs)+'\n'
    writelines+='module load '+modules+'\n\n'

    if short_recursive == 'True':
        time_lim = str(3600 * time - 100)
        # short_recursive command runs timer script before and after JDFT to check if walltime is hit
        writelines+='timeout 10 python '+os.path.join(jdftx_python_dir,'timer.py')+' > timer'+'\n'
        writelines+='timeout '+time_lim+' ' + 'python '+script+' > '+out+'\n'
        writelines+='timeout 10 python '+os.path.join(jdftx_python_dir,'timer.py')+' > timer'+'\n'

    else:
        writelines+='python '+script+' > '+out+'\n'
    writelines+='exit 0'+'\n'

    with open('submit.sh','w') as f:
        f.write(writelines)


def recursive_restart():
    with open('inputs','r') as f:
        inputs = f.read()
    lines = []
    for line in inputs.split('\n'):
        if 'restart' in line:
            if line == 'restart True':
                return
            lines.append('restart True')
        else:
            lines.append(line)
    with open('inputs','w') as f:
        f.write('\n'.join(lines))
    print('Updated inputs file for -r tag.')
    if os.path.exists('CONTCAR'):
        return
        #subprocess.call('cp CONTCAR CONTCAR_backup', shell=True)
    subprocess.call('cp POSCAR CONTCAR', shell=True)
    print('Copied POSCAR to CONTCAR')

if __name__ == '__main__':
    
    jdftx_python_dir = os.environ['JDFTx_manager_home']
    
    script = opj(jdftx_python_dir,'run_JDFTx.py')

    parser = argparse.ArgumentParser()
    parser.add_argument('-n', '--nodes', help='Number of nodes',
                        type=int, default=1)
    parser.add_argument('-c', '--cores', help='Number of cores',
                        type=int, default=os.environ['CORES_PER_NODE'])
    parser.add_argument('-t', '--time', help='Time limit',
                        type=int, default=int(os.environ['JDFTx_Default_Time']))
    parser.add_argument('-o', '--outfile', help='Outfile name',
                        type=str, required=True)
    parser.add_argument('-a', '--allocation', help='Allocation',
                        type=str,default='environ')
    parser.add_argument('-q', '--qos', help='Priority / QOS (e.g. high)',
                        type=str,default='standard')
    parser.add_argument('-r', '--recursive', help='If True, recursively runs job until complete',
                        type=str, default='False')
    parser.add_argument('-p', '--partition', help='Partition for Bridges2 (RM, RM-shared)',
                        type=str,default='RM-shared')
    parser.add_argument('-m', '--processes', help='Number of JDFT processes, should be <= nstates (see irr. kpoints). '+
                        'Total cores / processes = threads per process (int for high eff.)', type=int, default=2)
    parser.add_argument('-g', '--gpu', help='If True, runs GPU install of JDFTx',
                        type=str, default='False')
    parser.add_argument('-test', '--test_queue', help='If True, runs job on test queue',
                        type=str, default='False')

    args = parser.parse_args()
    #print(args.gpu, type(args.gpu))
    
    if args.recursive == 'True':
        recursive_restart()

    # new safety check line
    outfile = args.outfile
    if args.outfile == 'out':
        outfile = 'python_out'

    # Multiple write options depending on computer
    if comp == 'Eagle' or comp == 'Summit':
        write(args.nodes,args.cores,args.time,args.outfile,args.allocation,args.qos,		
              script, args.recursive, args.processes, args.gpu, args.test_queue)
    elif comp == 'Bridges2':
        write_bridges(args.nodes,args.cores,args.time,outfile,args.partition,args.qos,
                      script, args.recursive, args.processes)
    
    #os.system('sbatch submit.sh')
    subprocess.call('sbatch submit.sh', shell=True)
    
