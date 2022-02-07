
# Generating test data

I'm using a script to template the input files and to organize the 
files for easier version control. Each calculation is a `.yml` file.

## Generate the input files from the `.yml` file.

```sh
python3 ../__main__.py si init run --template_dir ../template_si --pseudo_dir ../pseudo
```

## Run the calculations.

### SCF

```sh
python3 ../__main__.py si scf run --nproc=8
```

### NSCF

```sh
python3 ../__main__.py si nscf run --nproc=8
```


### PW2WAN

```sh
python3 ../__main__.py si pw2wan run --nproc=8
```


### WANNIER90

```sh
python3 ../__main__.py si wannier90 run --nproc=8
```


