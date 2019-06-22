Quickstart
==========

bedshape is composed of two subcommands, *profile* and *alias*.

Using the profile subcommand
----------------------------

profile is the subcommand handling the generation of mutational and reactivity
profiles, plots; and *shape* files.

.. code-block:: bash

    $ bedshape profile --reference /path/to/reference.fa \
            --modified modified.bam \
            --region chr11:65505800-65505900

Like shapemapper2, bedshape requires only a modified sample to work with.
However, you may also pass a unmodified/untreated sample, and a denatured sample
as well.

.. code-block:: bash

    $ bedshape profile --reference /path/to/reference.fa \
            --modified modified.bam \
            --unmodified unmodified.bam \
            --denatured denatured.bam \
            --region chr11:65,505,800-65,505-900

The optional argument for the unmodified/untreated sample can be spelt as either
``--unmodified`` or ``--untreated``. Note that commas are allowed in the region,
as long as they come after the colon (``:``).

In lieu of a single region, bedshape can also take a BED file to produce outputs
for all regions specified therein.

.. code-block:: bash

    $ bedshape profile --reference /path/to/reference.fa \
            --modified modified.bam \
            --unmodified unmodified.bam \
            --bed interesting-genes.bed

Using the alias subcommand
--------------------------

You may often find yourself working with the same set of aligned reads, which
reside in some far-away directory with a very long name. To avoid having to
specify the location of these sets of files every time, bedshape supports the
*aliasing* of these files.

.. code-block:: bash

    $ bedshape alias add october-exp \
            --reference /path/to/reference.fa \
            --modified /data/october/modified.sam \
            --unmodified /data/october/unmodified.sam

Then, you may pass this 'alias' in place of ``--modified``, ``-unmodified``, and
``-denatured`` arguments in the profile subcommand.

.. code-block:: bash

    $ bedshape profile --alias october-exp \
            --bed interesting-genes.bed
