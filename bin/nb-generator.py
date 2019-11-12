#!/usr/bin/env python3
#
# nb-generator
#

import click

from pathlib import Path


def arg_parser():
    from argparse import ArgumentParser

    parser = ArgumentParser()
    return parser


@click.group()
def main():
    pass


@main.command()
@click.argument("filename")
@click.option("-o", "--output", default=None)
def mrc_view(filename, output):

    from nbformat import v4 as nbf, write as nbf_write

    from ROOT import TFile
    tfile = TFile.Open(filename)
    if not tfile:
        return 1

    nb = nbf.new_notebook()
    cells = nb['cells']
    cells.extend([
        nbf.new_markdown_cell(source=f"# MRC View {filename!r}"),
        nbf.new_code_cell('\n'.join([
            "import ROOT",
            "from ROOT import TFile, TCanvas"])),
        nbf.new_code_cell('\n'.join([
            "tfile = TFile.Open(%r)" % filename,
            "assert tfile"])),
        nbf.new_code_cell('\n'.join([
            "tdir = tfile.Get(%r)" % filename,
            "assert tdir"])),
            ""])),
    ])

    output_path = Path(output)
    with output_path.open("w") as f:
        nbf_write(nb, f)


if __name__ == "__main__":
    exit(main())
