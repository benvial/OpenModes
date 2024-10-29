
PROJECT_NAME := "openmodes"

BRANCH := "$(git branch --show-current)"

PROJECT_DIR := "$(realpath $PWD)"

VERSION := """$(python3 -c "import toml; print(toml.load('pyproject.toml')['project']['version'])")"""

version:
    @echo {{VERSION}}


doc:
    cd doc && make html

show:
    firefox doc/_build/html/index.html




all-fortran: set bld test-fortran

wrp:
    python -m numpy.f2py src/common.f90 src/rwg.f90 -m core -h core.pyf --lower \
    --f2cmap src/.f2py_f2cmap --backend meson \
    --build-dir build --overwrite-signature only:  \
    set_threads get_threads face_integrals_hanninen z_efie_faces_self \
    z_efie_faces_mutual arcioni_singular z_mfie_faces_self z_mfie_faces_mutual \
    face_integrals_yla_oijala :
    python -m numpy.f2py src/dunavant.f90 -m dunavant -h dunavant.pyf --lower \
    --f2cmap src/.f2py_f2cmap --backend meson \
    --build-dir build --overwrite-signature only: dunavant_order_num dunavant_rule :
    
set:
    meson setup --wipe builddir

bld:
    meson compile -C builddir

test-fortran:
    cd builddir && python -c "import openmodes"
    cd builddir && python -c "import openmodes.core"
    cd builddir && python -c "import openmodes.dunavant"

lspy:
    python -c 'import os; os.chdir("openmodes"); print([os.path.join("openmodes",f) for f in os.listdir() if os.path.isfile(f)])'
    python -c 'import os; os.chdir("openmodes/operator"); print([os.path.join("openmodes/operator",f) for f in os.listdir() if os.path.isfile(f)])'
    python -c 'import os; os.chdir("openmodes/external"); print([os.path.join("openmodes/external",f) for f in os.listdir() if os.path.isfile(f)])'
    python -c 'import os; os.chdir("openmodes/templates"); print([os.path.join("openmodes/templates",f) for f in os.listdir() if os.path.isfile(f)])'
    python -c 'import os; os.chdir("openmodes/static"); print([os.path.join("openmodes/static",f) for f in os.listdir() if os.path.isfile(f)])'
    python -c 'import os; os.chdir("openmodes/mesh"); print([os.path.join("openmodes/mesh",f) for f in os.listdir() if os.path.isfile(f)])'
    python -c 'import os; os.chdir("openmodes/geometry"); print([os.path.join("openmodes/geometry",f) for f in os.listdir() if os.path.isfile(f)])'



clean:
    rm -rf builddir build .coverage *.egg-info .pytest_cache htmlcov .ruff_cache

# Push to github
gh:
    @git add -A
    @read -p "Enter commit message: " MSG; \
    git commit -a -m "$MSG"
    @git push origin {{BRANCH}}


# Clean, reformat and push to github
save: gh

