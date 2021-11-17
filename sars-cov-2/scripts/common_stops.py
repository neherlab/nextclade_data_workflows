import pandas as pd
import click

@click.command()
@click.option('--number', default=20)
@click.option('--input-file', default="pre-processed/stops.tsv", help='Input file')
def format(number, input_file):
    df = pd.read_csv(input_file, sep="\t")

    output = []

    for row in df[1:number].sort_values(by=['stops']).itertuples():
        stop = row.stops.split(":")
        gene = stop[0]
        # beware: subtract 1 because in qc.json it's 0 indexed as opposed to nextclade output
        position = int(stop[1][1:-1])-1
        # {"geneName": "ORF7a", "codon": 61},

        output.append(f'{{ "geneName": "{gene}", "codon": {position} }},')
    click.echo("\n".join(sorted(output)))

if __name__ == '__main__':
    format()