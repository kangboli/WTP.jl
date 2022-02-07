import argparse
import subprocess
import re
from os import getcwd, path, mkdir, symlink, remove, listdir, environ, pathsep, chdir, cpu_count, removedirs, unlink
from shutil import copy, move, rmtree
from mako.template import Template
from yaml import safe_load, safe_dump
from time import sleep

wannier90_binary = '/global/homes/k/kangboli/.local/bin/wannier90.x'
parser = argparse.ArgumentParser(description='Run a quantum chem calculation')

parser.add_argument('prefix', type=str,
                    help='The prefix of the system. Same as the prefix in quantum espresso input files')
parser.add_argument('calculation', type=str,
                    help='The type of calculation. can be'
                         'init; scf; nscf; pp; pw2wan; wannier90')
parser.add_argument('op', type=str,
        help='The operation to be performed regarding the calculation. can be: '
                         'run; clean')
parser.add_argument('--dir', dest='dir',
                    default=getcwd(),
                    help='The directory that has the folder input and output. Use the current dir by default.'
                         'This is only for running code locally.'
                         'For running code on the cluster, the input dir and the output dir are different.')
parser.add_argument('--nproc', dest='nproc', type=int, default=cpu_count(),
        help='Local only: The number of processes to use.')
parser.add_argument('--cluster', action='store_true',
                    help='Whether or not the code is running on a cluster.')
parser.add_argument('--wall_time', dest='wall_time', default='00:30:00',
                    help='Cluster only: The wall time for the calculation.')
parser.add_argument('--qos', dest='qos', default='regular',
        help='Cluster only: The options are regular (default), flex, debug.')
parser.add_argument('--nodes', dest='nodes', type=int, default=1,
        help='Cluster only: The number of nodes to use.')
parser.add_argument('--tasks_per_node', dest='tpn', type=int, default=32,
        help='Cluster only: The number of nodes to use.')
parser.add_argument('--template_dir', dest='template_dir', type=str,
        help='Initialize the input files from a template folder. They must have the file extention:\n'
        '.scf .nscf .pw2wan and .win.')
parser.add_argument('--pseudo_dir', dest='pseudo_dir', type=str,
        help='The directory with the pseudo potentials.')

args = parser.parse_args()
prefix = args.prefix
prefix_dir = path.realpath(args.dir)
chdir(prefix_dir)
on_cluster = args.cluster
wall_time = args.wall_time
qos = args.qos
nproc = args.nproc
nodes = args.nodes
tpn = args.tpn
template_dir = args.template_dir
pseudo_dir = args.pseudo_dir

if args.calculation == 'init' and not pseudo_dir and args.op != 'clean':
    print('--pseudo_dir is required at init.')
    exit()

if tpn < 32 and nodes > 1:
    print('Don\'t use more than one node unless there are 32 tasks per node.')
    exit()

if args.calculation == 'pp' and (nodes > 1 or tpn > 1) and on_cluster:
    print(f'wannier90.x -pp should not run in parallel: {nodes} nodes and {tpn} tasks per nodes have been specified')
    exit()

in_dir = path.join(prefix_dir, 'input')
scf_in = path.join(in_dir, f'{prefix}.scf')
nscf_in = path.join(in_dir, f'{prefix}.nscf')
pp_in = path.join(in_dir, f'{prefix}.win')
wannier_in = path.join(in_dir, f'{prefix}.win')
pw2wan_in = path.join(in_dir, f'{prefix}.pw2wan')
yaml_file = path.join(prefix_dir, f'{prefix}.yaml')

out_dir = path.join(prefix_dir, 'output')
scf_out = path.join(prefix_dir, 'output', 'scf')
nscf_out = path.join(prefix_dir, 'output', 'nscf')
wfc_out = path.join(prefix_dir, 'output', 'wfc')
wannier90_out = path.join(prefix_dir, 'output', 'win')
pw2wan_out = path.join(prefix_dir, 'output', 'pw2wan')
scdm_out = path.join(prefix_dir, 'output', 'scdm')

scf_out_files = [path.join(scf_out, f'{prefix}{postfix}') for postfix in ['_scf.out', '.save', '.xml']]
nscf_out_files = [path.join(nscf_out, f'{prefix}{postfix}') for postfix in ['_nscf.out', '.save', 'xml']]
pp_out_files = [path.join(wannier90_out, f'{prefix}{postfix}') for postfix in ['.nnkp', '.wout']]
pw2wan_out_files = [path.join(pw2wan_out, f'{prefix}{postfix}') for postfix in
                    ['_pw2wan.out', '.amn', '.eig', '.mmn']]
wannier90_out_files = [path.join(wannier90_out, f'{prefix}.wout')]

# time_switch = '-f'
# time_flag = '\"CPU Time: Kernel -> %Ss  User -> %Us; Max Memory: %M KB; IO: In -> %I Out -> %O\"'
# mpi_flags = '--use-hwthread-cpus'
time_switch = ''
time_flag = ''
mpi_flags = ''

parallel = ['mpirun', mpi_flags, '-np', str(nproc)]

yml = open(yaml_file, 'r')
config = safe_load(yml)
yml.close()

# For Clusters

def remove_file(filename):
    try:
        if path.isfile(filename):
            remove(filename)
        else:
            rmtree(filename)
    except FileNotFoundError as e:
        print(f"{filename} ignored since it does not exist.")

def move_overwrite(src: str, dst):
    if path.isdir(dst):
        dst = path.join(dst, src.split('/')[-1])
    if path.exists(dst):
        print(f'Existing {dst} will be overwritten.')
        remove_file(dst)
    # print(f'Moving {src} to {dst}')
    move(src, dst)

def get_dependency():
    if args.calculation == 'scf':
        return ''
    elif args.calculation == 'nscf':
        return f'--dependency=afterany:{status_query("scf")["id"]}'
    elif args.calculation == 'pp':
        return f'--dependency=afterany:{status_query("nscf")["id"]}'
    elif args.calculation == 'pw2wan':
        return f'--dependency=afterany:{status_query("pp")["id"]}'
    elif args.calculation == 'wannier90':
        return f'--dependency=afterany:{status_query("pw2wan")["id"]}'

def submit_job(binary: str, wall_time: str, input_file=None, output_file=None, tasks_per_node=32):

    while get_dependency() == 'nil':
        sleep(1)

    job_script_wannier = f'#!/bin/bash -l\n' \
                     f'#SBATCH -q {qos}\n' \
                     f'#SBATCH --time={wall_time}\n' \
                     f'#SBATCH --nodes={nodes}\n' \
                     f'#SBATCH --constraint=haswell\n' \
                     f'#SBATCH --tasks-per-node={tasks_per_node}\n' \
                     f'#SBATCH --cpus-per-task=2\n' \
                     f'module purge\n' \
                     f'module load PrgEnv-intel/6.0.5\n' \
                     f'module load impi/2020\n' \
                     f'srun {get_dependency()} --cpu-bind=cores {binary}\n'

# -n {nodes*tasks_per_node} -N {nodes} -c 1

    job_script_mpi = f'#!/bin/bash -l\n' \
                     f'#SBATCH -q {qos}\n' \
                     f'#SBATCH --time={wall_time}\n' \
                     f'#SBATCH --nodes={nodes}\n' \
                     f'#SBATCH --constraint=haswell\n' \
                     f'#SBATCH --tasks-per-node={tasks_per_node}\n' \
                     f'#SBATCH -A m3799\n' \
                     f'export OMP_PROC_BIND=true\n' \
                     f'export OMP_PLACES=cores\n' \
                     f'export OMP_NUM_THREADS=1\n' \
                     f'module purge\n' \
                     f'module load PrgEnv-intel/6.0.5\n' \
                     f'module load impi/2020\n' \
                     f'srun {get_dependency()} -n {nodes*tasks_per_node} -N {nodes} -c 1 {binary} < {input_file} > {output_file}\n'

    submission_file = f'{prefix_dir}/submission_{args.calculation}.sh'
    with open(submission_file, 'w+') as submission:
        if 'wannier90' in binary and 'pw2wannier90' not in binary:
            submission.write(job_script_wannier)
        else:
            submission.write(job_script_mpi)

    cmd = ['sbatch', submission_file]
    result = subprocess.run(cmd, capture_output=True, encoding='UTF-8')
    m = re.match(r'Submitted\s+batch\s+job\s+(\d+)', result.stdout)
    job_id = m.group(1)
    print(result.stdout)
    print(result.stderr)
    status_update(args.calculation, job_id, 'running')
    # remove_file(f'{prefix_dir}/submission.sh')
    # Wait for completion.
    
    queue = job_id
    while job_id in queue:
        result = subprocess.run(['squeue', '-u', 'kangboli'], capture_output=True, encoding='UTF-8')
        queue = result.stdout
        sleep(3)
    
    status_update(args.calculation, job_id, 'done')
    print('Job completed')
    

sys_path = environ['PATH'].split(pathsep)

def clean_scf():
    remove_file(f'{prefix}.save')
    remove_file(f'{prefix}.xml')
    for file in scf_out_files:
        remove_file(file)

def clean_nscf():
    remove_file(f'{prefix}.save')
    remove_file(f'{prefix}.xml')
    for file in nscf_out_files:
        remove_file(file)

def clean_pp():
    for file in pp_out_files:
        remove_file(file)
    remove_file(f'{prefix}.nnkp')

def clean_pw2wan():
    if path.exists('./unk'):
        unlink('./unk')
    for file in pw2wan_out_files:
        remove_file(file)


def clean_wannier90():
    for file in wannier90_out_files:
        remove_file(file)


def clean_init():
    if path.exists('./pseudo'):
        unlink(f'./pseudo')
    if path.exists('./unk'):
        unlink('./unk')
    for file in listdir(prefix_dir):
        if not (str.endswith(file, 'yml') or str.endswith(file, 'yaml') or
                str.endswith(file, 'org') or str.endswith(file, 'md')):
            remove_file(file)
        if path.exists('status.yaml'):
            remove_file('status.yaml')

if args.op == 'clean':
    if args.calculation == 'scf':
        clean_scf()
    if args.calculation == 'nscf':
        clean_nscf()
    if args.calculation == 'pp':
        clean_pp()
    if args.calculation == 'pw2wan':
        clean_pw2wan()
    if args.calculation == 'wannier90':
        clean_wannier90()
    if args.calculation == 'init':
        clean_init()
    exit()

def find_binary(name):

    for p in sys_path:
        if path.exists(p) and name in listdir(p):
            return path.join(p, name)
    else:
        print(f'{name} is not on system PATH')
        exit(1)


def create_dirs():
    for d in [in_dir, out_dir, wfc_out, scf_out, nscf_out, wannier90_out, pw2wan_out, scdm_out, '{}/unk'.format(pw2wan_out)]:
        if not path.exists(d) or not path.isdir(d):
            mkdir(d)


def check_dirs():
    for d in [in_dir, out_dir, wfc_out, scf_out, nscf_out, wannier90_out, pw2wan_out, scdm_out, '{}/unk'.format(pw2wan_out)]:
        if not path.exists(d) or not path.isdir(d):
            print(f'{d} does not exist. Please init first')
            exit()


def render_input_file(filename):
    template = Template(filename=filename)
    extension = filename.split('.')[-1]
    with open(f'{in_dir}/{prefix}.{extension}', 'w+') as out_file:
        out_file.write(template.render(**config))
    
def status_update(calculation, job_id, status):
    with open(f'{prefix_dir}/status.yaml', 'r') as yml:
        cal_status = safe_load(yml)
    cal_status[calculation] = {'id': job_id, 'status': status}
    with open(f'{prefix_dir}/status.yaml', 'w+') as yml:
        yml.write(safe_dump(cal_status))

def status_query(calculation):
    with open(f'{prefix_dir}/status.yaml') as yml:
        cal_status = safe_load(yml)
    return cal_status[calculation]


if args.calculation == 'init':
    create_dirs()
    if len(template_dir) > 0:
        for file in listdir(template_dir):
            render_input_file(path.realpath(path.join(template_dir, file)))
    if on_cluster:
        cal_status = {'scf': {'id': 'nil', 'status': '' }, 
                'nscf': {'id': 'nil', 'status': '' }, 
                'pp': {'id': 'nil', 'status': '' }, 
                'pw2wan': {'id': 'nil', 'status': '' },  
                'wannier90': {'id': 'nil', 'status': '' }}
        with open(f'{prefix_dir}/status.yaml', 'w+') as yml:
            yml.write(safe_dump(cal_status))

    try:
        symlink(pseudo_dir, f'{prefix_dir}/pseudo')
    except FileExistsError as e:
        pass
        
else:
    check_dirs()


def run_cmd_with_input(cmd, input_file, output_file):
    print('Running:')
    print(f'{" ".join(cmd)} < {input_file} > {output_file}')
    cmd = [c for c in cmd if len(c) > 0]
    with open('cmd.sh', 'a') as cmd_out:
        cmd_out.write(' '.join(cmd) + '\n')
    if not path.exists(input_file):
        print(f'{input_file} is required.')
    with open(input_file, 'r') as fin:
        result = subprocess.run(cmd, input=fin.read(), capture_output=True, encoding='UTF-8')
    with open(output_file, 'w') as out:
        out.write(result.stdout)
    if result.stderr != '':
        print(result.stderr)


if args.calculation not in ['init', 'scf', 'nscf', 'pp', 'pw2wan', 'wannier90']:
    print(f'{args.calculation} is not supported.') 
    exit()

# Need to check that hte number of processors match between scf and nscf.

if args.calculation == 'scf' or args.calculation == 'nscf':
    pwx = [find_binary('pw.x'), '-nk', str(nproc), '-nb', str(1)]
    cmd = ['/usr/bin/time', time_switch, time_flag] + parallel + pwx
    input_file = scf_in if args.calculation == 'scf' else nscf_in
    calc_out = scf_out if args.calculation == 'scf' else nscf_out
    if on_cluster:
        submit_job(' '.join(cmd), wall_time,
                   input_file, path.join(calc_out, f'{prefix}_{args.calculation}.out'), tasks_per_node=tpn)
    else:
        run_cmd_with_input(
            cmd, input_file, path.join(calc_out, f'{prefix}_{args.calculation}.out'))
    if args.calculation == 'nscf':
        for file in listdir(prefix_dir):
            if 'wfc' in file:
                move_overwrite(file, wfc_out)
    try:
        symlink(path.join(prefix_dir, f'{prefix}.save'), path.join(calc_out, f'{prefix}.save'))
        symlink(path.join(prefix_dir, f'{prefix}.xml'), path.join(calc_out, f'{prefix}.xml'))
    except FileExistsError as e:
        pass


if args.calculation == 'pp':
    if path.exists(wannier90_binary):
        cmd = [wannier90_binary, '-pp', prefix]
    else:
        cmd = [find_binary('wannier90.x'), '-pp', prefix]
        
    copy(pp_in, f'{prefix}.win')

    with open('cmd.sh', 'a') as cmd_out:
        cmd_out.write(' '.join(cmd) + '\n')

    if on_cluster:
        submit_job(' '.join(cmd), wall_time, tasks_per_node=1)
    else:
        subprocess.run(cmd)
    remove_file(f'{prefix}.win')
    copy(f'{prefix}.nnkp', wannier90_out)

    copy(f'{prefix}.wout', wannier90_out)
    remove_file(f'{prefix}.wout')


if args.calculation == 'pw2wan':
    cmd = ['/usr/bin/time', time_switch, time_flag] + parallel + [find_binary('pw2wannier90.x')]
    if on_cluster:
        submit_job(find_binary('pw2wannier90.x'), wall_time,
                   input_file = pw2wan_in, output_file = path.join(pw2wan_out, f'{prefix}_{args.calculation}.out'), tasks_per_node=tpn)
    else:
        run_cmd_with_input(cmd, pw2wan_in, path.join(pw2wan_out, f'{prefix}_{args.calculation}.out'))

    for postfix in ['amn', 'eig', 'mmn']:
        file = path.join(prefix_dir, f'{prefix}.{postfix}')
        move_overwrite(file, pw2wan_out)

    for file in listdir('.'):
        if file.startswith('UNK'):
            move_overwrite(file, path.join(pw2wan_out, 'unk'))
    if path.exists('./unk'):
        unlink('./unk')
    symlink(path.join(pw2wan_out, 'unk'), './unk')

    remove_file('input_tmp.in')

if args.calculation == 'wannier90':
    copy(wannier_in, f'{prefix}.win')
    try:
        move_overwrite(f'{pw2wan_out}/{prefix}.amn', f'{prefix}.amn')
        move_overwrite(f'{pw2wan_out}/{prefix}.eig', f'{prefix}.eig')
        move_overwrite(f'{pw2wan_out}/{prefix}.mmn', f'{prefix}.mmn')
    except FileNotFoundError as e:
        pass

    if on_cluster:
        if path.exists(wannier90_binary):
            submit_job(wannier90_binary + ' ' + prefix, wall_time, tasks_per_node=tpn)
        else:
            submit_job(find_binary('wannier90.x') + ' ' + prefix, wall_time, tasks_per_node=tpn)
    else:
        cmd = ['/usr/bin/time', time_switch, time_flag, find_binary('wannier90.x'), prefix]
        cmd = [c for c in cmd if len(c) > 0]
        with open('cmd.sh', 'a') as cmd_out:
            cmd_out.write(' '.join(cmd) + '\n')
        subprocess.run(cmd, encoding='UTF-8')
    move_overwrite(f'{prefix}.wout', f'{wannier90_out}/{prefix}.wout')
    move_overwrite(f'{prefix}_centres.xyz', f'{wannier90_out}/{prefix}_centres.xyz')
    move_overwrite(f'{prefix}_u.mat', f'{wannier90_out}/{prefix}_u.mat')

    move_overwrite(f'{prefix}.amn', f'{pw2wan_out}/{prefix}.amn')
    move_overwrite(f'{prefix}.eig', f'{pw2wan_out}/{prefix}.eig')
    move_overwrite(f'{prefix}.mmn', f'{pw2wan_out}/{prefix}.mmn')

    remove(f'{prefix}.win')

