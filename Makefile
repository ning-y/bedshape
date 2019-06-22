BINARIES = src/shapemapper_mutation_parser \
        src/shapemapper_mutation_counter \
        src/make_reactivity_profiles.py \
        src/normalize_profiles.py \
        src/render_figures.py \
        src/tab_to_shape.py

all: $(BINARIES)
	pipenv install

docs: sphinx src/cli.py
	pipenv run sphinx-build -M html sphinx sphinx-docs
	rm -rf docs && mkdir docs && mv sphinx-docs/html/* docs/
	echo "bedshape.ningyuan.io" > docs/CNAME
	touch docs/.nojekyll

src/shapemapper_mutation_parser: shapemapper2
	bash shapemapper2/internals/install/build_binaries.sh
	cp shapemapper2/internals/bin/shapemapper_mutation_parser src/
src/shapemapper_mutation_counter: src/shapemapper_mutation_parser
	cp shapemapper2/internals/bin/shapemapper_mutation_counter src/

src/make_reactivity_profiles.py: shapemapper2
	cp shapemapper2/internals/bin/make_reactivity_profiles.py src/
src/normalize_profiles.py: shapemapper2
	cp shapemapper2/internals/bin/normalize_profiles.py src/
src/render_figures.py: shapemapper2
	cp shapemapper2/internals/bin/render_figures.py src/
src/tab_to_shape.py: shapemapper2
	cp shapemapper2/internals/bin/tab_to_shape.py src/

shapemapper2:
	git submodule init
	git submodule update
