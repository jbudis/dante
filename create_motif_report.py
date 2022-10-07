import os
import argparse
import textwrap

import numpy as np
import plotly.graph_objects as go
from plotly.io import to_json
from bs4 import BeautifulSoup

motif_summary = """
<h2 id="summary">Summary table</h2>
<table class="tg" id="tg">
    <thead>
        <tr>
            <th class="tg-s6z2" rowspan="2">Sample</th>
            <th class="tg-s6z2" colspan="2">Allele 1</th>
            <th class="tg-s6z2" colspan="2">Allele 2</th>
            <th class="tg-s6z2" rowspan="2">Overall<br>confidence</th>
            <th class="tg-s6z2" colspan="2">Reads</th>
        </tr>
        <tr>
            <td class="tg-s6z2">prediction</td>
            <td class="tg-s6z2">confidence</td>
            <td class="tg-s6z2">prediction</td>
            <td class="tg-s6z2">confidence</td>
            <td class="tg-s6z2">full</td>
            <td class="tg-s6z2">partial</td>
        </tr>
    </thead>
    <tbody>
        {table}
    </tbody>
</table>
<script>
    $(document).ready( function () {{
    $('#tg').DataTable();
}} );
</script>
"""

row_string = """  <tr>
    <td class="tg-s6z2">{name}</td>
    <td class="tg-s6z2">{allele1}</td>
    <td class="tg-s6z2">{allele1_conf}%</td>
    <td class="tg-s6z2">{allele2}</td>
    <td class="tg-s6z2">{allele2_conf}%</td>
    <td class="tg-s6z2">{motif_conf}%</td>
    <td class="tg-s6z2">{reads_blue}</td>
    <td class="tg-s6z2">{reads_grey}</td>
  </tr>
"""

heatmap_string = """
<div id="motif-heatmap"></div>
<script>
    Plotly.newPlot('motif-heatmap', {heatmap}, {{}});
</script>
"""


def load_arguments():
    """
    Loads and parses arguments.
    :return: input_dir - path to directory with Dante reports
             output_dir - path to output dir
             args - parsed arguments
    """
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description=textwrap.dedent("""Python program to collect Dante report files and 
                                     create 'reversed' report for each motif, to make the comparison of multiple samples 
                                     easier"""))

    # add arguments
    parser.add_argument('input_dir', help='Path to directory with Dante reports')
    parser.add_argument('output_dir', help='Path to directory where the output will be stored', nargs='?',
                        default='example/motif_report')

    args = parser.parse_args()

    input_dir = args.input_dir
    output_dir = args.output_dir

    return input_dir, output_dir, args


def custom_format(template, **kwargs):
    """
    Custom format of strings for only those that we provide
    :param template: str - string to format
    :param kwargs: dict - dictionary of strings to replace
    :return: str - formatted string
    """
    for k, v in kwargs.items():
        template = template.replace('{%s}' % k, v)

    return template


def generate_row(motif_name, a1, c1, a2, c2, c, reads_blue, reads_grey):
    return row_string.format(name=motif_name, allele1=a1, allele1_conf=c1, allele2=a2, allele2_conf=c2,
                             motif_conf=c, reads_blue=reads_blue, reads_grey=reads_grey)


# Convert extracted numbers from table to ints and set background and expanded alleles to 0
def parse_alleles(num):
    if num == 'B' or num == 'E':
        return 0
    else:
        return int(num)


def parse_label(num):
    if num == 0:
        return ''
    else:
        return str(int(num))


def generate_motif_report(path, key, samples, fig):
    """
    Generate report file for one motif
    :param path: str - path to output dir
    :param key: str - motif name
    :param samples: list - list of samples of selected motif
    :param fig: str - heatmap object
    """
    template = open('report/motif_report.html', 'r').read()
    rows = []
    key = key.replace('/', '-')

    for sample in samples:
        rows.append(generate_row(sample[0], sample[1], sample[2], sample[3], sample[4], sample[5], sample[6], sample[7]))

    with open('%s/report_%s.html' % (path, key), 'w') as f:
        table = motif_summary.format(table='\n'.join(rows))
        heatmap = heatmap_string.format(heatmap=to_json(fig))
        f.write(custom_format(template, motif=key.split('_')[0], seq=key.split('_')[1], motifs_content=table, motif_heatmap=heatmap))


def create_reports(input_dir, output_dir):
    os.makedirs(output_dir, exist_ok=True)
    paths = []

    # find all report files in root directory
    for root, dirs, files in os.walk(input_dir):
        for file in files:
            if file.endswith('.html'):
                paths.append(os.path.join(root, file))

    motifs = {}

    for path in paths:
        file = BeautifulSoup(open(path, 'r'), 'html.parser')
        fname = path.split('/')[-1].split('.')[0]

        # find table of class 'tg' and extract all rows from it
        for row in file.find(class_='tg').find_all('tr'):
            columns = row.find_all('td')

            # remove head rows
            if columns == [] or columns[0].text.strip() == 'prediction':
                continue

            if columns:
                name = columns[0].text.strip() + '_' + columns[1].text.strip()
                doc = [fname, columns[2].text.strip(), columns[3].text.strip().replace('%', ''),
                       columns[4].text.strip(), columns[5].text.strip().replace('%', ''),
                       columns[6].text.strip().replace('%', ''), columns[7].text.strip(), columns[8].text.strip()]

                if name not in motifs:
                    motifs[name] = [doc]
                else:
                    motifs[name].append(doc)

    # create heatmap of alleles
    for _key in motifs.keys():
        a1 = [parse_alleles(row[1]) for row in motifs[_key]]
        a2 = [parse_alleles(row[3]) for row in motifs[_key]]

        arr = np.zeros((max(a1) + 1, max(a2) + 1))

        for i in range(len(a1)):
            arr[a1[i], a2[i]] += 1

        text = [[parse_label(arr[i, j]) for j in range(arr.shape[1])] for i in range(arr.shape[0])]

        arr = arr / np.max(arr)

        fig = go.Figure(go.Heatmap(z=list(arr), text=text, texttemplate="%{text}", textfont={"size": 10},
                                   colorscale='Hot_r'))
        fig.update_layout(width=750, height=750, template='simple_white')
        fig.update_yaxes(title_text="Allele 1")
        fig.update_xaxes(title_text="Allele 2")

        # generate motif
        generate_motif_report(output_dir, _key, motifs[_key], fig)


if __name__ == '__main__':
    input_d, output_d, args = load_arguments()
    create_reports(input_d, output_d)
