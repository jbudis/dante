import os
import argparse
import re
import textwrap
from itertools import chain

import numpy as np
import plotly.graph_objects as go
from plotly.io import to_json
from bs4 import BeautifulSoup
from natsort import natsorted

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
    $('#tg').DataTable( {{
        'order': [],
        dom: 'Bfrtip',
        buttons: [
            'copy', 'csv', 'excel', 'pdf', 'print'
        ]
    }} );
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

summary_plot_string = """
<h2 id="map">Allele heatmap</h2>
<p><b>Number of background results: {background}</b></p>
<div id="motif-heatmap"></div>
<script>
    Plotly.newPlot('motif-heatmap', {heatmap}, {{}});
</script>
<h2 id="hist">Allele histogram</h2>
<div id="motif-hist"></div>
<script>
    Plotly.newPlot('motif-hist', {histogram}, {{}});
</script>
"""

motif_plot_string = """<h3>{sample_name}</h3>
<table>
    <tr>
        <td colspan="2">
            <div class="hist pic100" id="hist-{sample_id}"></div>
            <script>
                Plotly.react('hist-{sample_id}', {hist_plot}, {{}});
            </script>
        </td>
        <td colspan="1">
            <div class="pcol pic100" id="pcol-{sample_id}"></div>
            <script>
                Plotly.react('pcol-{sample_id}', {pcol_plot}, {{}});
            </script>
        </td>
    </tr>
</table>
"""

align_string = """
  <details>
    <summary>{display_text}</summary>
    <div id="align-{sample_id}" class="align">press "Run with JS"</div>
    <script>
        var fasta = {fasta};
        var seqs = msa.io.fasta.parse(fasta);
        var opts = {{
            el: document.getElementById("align-{sample_id}"),
            vis: {{
                conserv: false,
                metaIdentity: true,
                overviewbox: true,
                seqlogo: true
            }},
            seqs: seqs,
            colorscheme: {{"scheme": "nucleotide"}},
            // smaller menu for JSBin
            menu: "small",
            bootstrapMenu: true
        }};
        var m = new msa.msa(opts);
        m.render()
    </script>
  </details>
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
    parser.add_argument('--report_every', type=int, help='Specify how often a progress message should be printed (default=5)',
                        default=5)
    parser.add_argument('-q', '--quiet', help="Don't print any progress messages", action='store_true')

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
        return num
    elif num == '---':
        return -1
    else:
        return int(num)


def parse_label(num):
    if num == 0:
        return ''
    else:
        return str(int(num))


def generate_motif_report(path, key, samples, plots, alignments, fig_heatmap, fig_hist, bgs):
    """
    Generate report file for one motif
    :param path: str - path to output dir
    :param key: str - motif name
    :param samples: list - list of samples of selected motif
    :param plots: list - list of plots of selected motif
    :param alignments: list - list of alignments of selected motif
    :param fig_heatmap: str - heatmap object
    :param fig_hist: str - histogram object
    :param bgs: int - number of background results
    """
    script_dir = os.path.dirname(os.path.abspath(__file__))
    template = open('%s/report/motif_report.html' % script_dir, 'r').read()
    rows = []
    key = key.replace('/', '-')

    for sample in samples:
        rows.append(generate_row(sample[0], sample[1], sample[2], sample[3],
                                 sample[4], sample[5], sample[6], sample[7]))

    motif_plots = list(chain.from_iterable(zip(plots, alignments)))

    with open('%s/report_%s.html' % (path, key), 'w') as f:
        table = motif_summary.format(table='\n'.join(rows))
        summary_plots = summary_plot_string.format(background=bgs, heatmap=to_json(fig_heatmap), histogram=to_json(fig_hist))
        f.write(custom_format(template, motif=key.split('_')[0], seq=key.split('_')[1],
                              motifs_content=table, summary_plots=summary_plots, motif_plots='\n'.join(motif_plots)))


def create_reports(input_dir, output_dir, arg_list):
    """
    Traverse input dir, collect all tables and plots from reports and generate motif reports
    :param input_dir: str - path to input dir
    :param output_dir: str - path to output dir
    :param arg_list: [any] - list of arguments
    """
    os.makedirs(output_dir, exist_ok=True)
    paths = []

    # find all report files in root directory
    for root, dirs, files in os.walk(input_dir):
        for file in files:
            if file == 'report.html':
                paths.append(os.path.join(root, file))

    motifs = {}
    plots = {}
    alignments = {}
    current_alignment = []
    prev_name = ''

    if len(paths) == 0:
        print("ERROR\tInput directory is empty")
        return

    cnt = 1
    print("INFO\tParsing reports")

    for path in natsorted(paths):
        if not arg_list.quiet and cnt % arg_list.report_every == 0:
            print("INFO\tParsing file\t%d/%d" % (cnt, len(paths)))
        cnt += 1

        file = BeautifulSoup(open(path, 'r'), 'html.parser')
        sample = file.find(id='sample_name').text.strip()

        # find table of class 'tg' and extract all rows from it
        for cl in file.find_all(class_='tg'):
            for row in cl.find_all('tr'):
                columns = row.find_all('td')

                # remove head rows
                if columns == [] or columns[0].text.strip() == 'prediction':
                    continue

                if len(columns) >= 9:
                    name = re.sub(r'[^\w_]', '', columns[0].text.strip()) + '_' + columns[1].text.strip()
                    if ',' in name:
                        break

                    doc = [sample, columns[2].text.strip(), columns[3].text.strip().replace('%', ''),
                           columns[4].text.strip(), columns[5].text.strip().replace('%', ''),
                           columns[6].text.strip().replace('%', ''), columns[7].text.strip(), columns[8].text.strip()]

                    if name not in motifs:
                        motifs[name] = [doc]
                    else:
                        motifs[name].append(doc)

        for cl in file.find_all(class_='plots'):
            name = cl.find(class_='hist')['class'][-1]
            hist = cl.find(class_='hist').find_next('script').text
            pcol = cl.find(class_='pcol')

            if pcol:
                pcol = pcol.find_next('script').text

            prev = cl.find_previous('p').find_all('u')

            if len(prev) > 1:
                continue

            hist_data_re = re.match(r"[\w\W]+let hist_data = ({.+});[\w\W]+", hist)
            if hist_data_re is None:
                continue
            hist_data = hist_data_re.group(1)
            pcol_data_re = re.match(r"[\w\W]+let pcol_data = ({.+});[\w\W]+", pcol)
            if pcol_data_re is None:
                continue
            pcol_data = pcol_data_re.group(1)

            name += '_' + prev[0].text
            temp = motif_plot_string.format(sample_name=sample, sample_id=name + '_' + sample,
                                            hist_plot=hist_data, pcol_plot=pcol_data)

            if name not in plots:
                plots[name] = [temp]
            else:
                plots[name].append(temp)

        for cl in file.find_all(class_='align'):
            name = cl['class'][-1]
            msa = cl.find_next('script').text
            prev = cl.find_previous('p').find_all('u')
            disp = cl.find_previous('summary').text

            if len(prev) > 1:
                continue

            msa_data = re.match(r"[\w\W]+let \w+_fasta = (`[\w\W]*`);[\w\W]+", msa).group(1)

            name += '_' + prev[0].text
            if prev_name == '':
                prev_name = name

            temp = align_string.format(sample_name=sample, sample_id=cl['id'] + sample,
                                       fasta=msa_data, display_text=disp)

            if prev_name == name:
                current_alignment.append(temp)
            else:
                if prev_name not in alignments:
                    alignments[prev_name] = ['\n'.join(current_alignment)]
                else:
                    alignments[prev_name].append('\n'.join(current_alignment))

                current_alignment = [temp]
                prev_name = name

        if len(current_alignment) != 0:
            if prev_name not in alignments:
                alignments[prev_name] = ['\n'.join(current_alignment)]
            else:
                alignments[prev_name].append('\n'.join(current_alignment))

        current_alignment = []
        prev_name = ''


    print("INFO\tGenerating motif reports")

    # create histogram of read counts
    for _key in motifs.keys():
        a1 = [parse_alleles(row[1]) for row in motifs[_key] if parse_alleles(row[1]) != -1]
        a2 = [parse_alleles(row[3]) for row in motifs[_key] if parse_alleles(row[3]) != -1]

        if len(a1) == 0 or len(a2) == 0:
            continue

        try:
            a1_max = max(x for x in a1 if isinstance(x, int))
        except ValueError:
            a1_max = 0

        try:
            a2_max = max(x for x in a2 if isinstance(x, int))
        except ValueError:
            a2_max = 0

        arr = np.zeros((a1_max + 2, a2_max + 2))

        max_count = max(a1_max, a2_max)
        hist_arr = [0 for _ in range(max_count + 2)]
        hist_color = ['#636EFA' for _ in range(max_count + 1)] + ['#EF553B']

        bgs = 0

        for i in range(len(a1)):
            if a2[i] == 'E':
                if a1[i] == 'E':
                    arr[a1_max + 1, a2_max + 1] += 1
                else:
                    arr[a1[i], a2_max + 1] += 1
                hist_arr[-1] += 1
            elif a2[i] == 'B':
                bgs += 1
            elif a2[i] == '---':
                pass
            else:
                arr[a1[i], a2[i]] += 1
                hist_arr[a1[i]] += 1
                hist_arr[a2[i]] += 1

        fig_histogram = go.Figure(data=[
            go.Bar(y=hist_arr, text=[parse_label(num) for num in hist_arr], name='Count histogram', marker_color=hist_color),
        ])
        #fig_histogram.add_vline(x=len(hist_arr) - 1.5, line_width=5, line_color='black', opacity=1)
        fig_histogram.update_xaxes(title_text="Prediction", tickmode='array',
                                   tickvals=np.concatenate([np.array(range(0, max_count + 1, 5)), [max_count + 1]]),
                                   ticktext=list(range(0, max_count + 1, 5)) + ['E(>%d)' % (max_count + 1)])
        fig_histogram.update_yaxes(title_text="Count")
        fig_histogram.update_traces(hovertemplate="<b>Prediction:\t%{x}</b><br />Count:\t%{y}<br />", textfont_size=7)
        fig_histogram.update_layout(width=1000, height=500, template='simple_white',
                                    barmode='stack', yaxis_fixedrange=True, hovermode='x')

        row_max, column_max = arr.shape
        row_sum = np.sum(arr, axis=1)
        row_count, column_count = 0, 0

        for element in row_sum:
            if element == 0:
                row_count += 1
            else:
                break

        column_sum = np.sum(arr, axis=0)

        for element in column_sum:
            if element == 0:
                column_count += 1
            else:
                break

        text = [[parse_label(arr[i, j]) for j in range(arr.shape[1])] for i in range(arr.shape[0])]

        arr = arr / np.max(arr)

        # create heatmap of alleles
        fig_heatmap = go.Figure(data=[
            go.Heatmap(z=list(arr), text=text, textfont={"size": 10}, colorscale='Hot_r',
                       hovertemplate="<b>Allele 1:\t%{y}<br />Allele 2:\t%{x}</b><br />Count:\t%{text}",
                       texttemplate="%{text}", name='Prediction heatmap')
        ])
        fig_heatmap.add_vline(x=a2_max + 0.5, line_width=5, line_color='black', opacity=1)
        fig_heatmap.add_hline(y=a1_max + 0.5, line_width=5, line_color='black', opacity=1)
        fig_heatmap.update_layout(width=750, height=750, template='simple_white',
                                  yaxis=dict(range=[row_count-1.5, row_max-0.5]),
                                  xaxis=dict(range=[column_count-1.5, column_max-0.5]))
        fig_heatmap.update_yaxes(title_text="Allele 1", tickmode='array',
                                 tickvals=np.concatenate([np.array(range(0, row_max - 1, 5)), [row_max - 1]]),
                                 ticktext=list(range(0, row_max - 1, 5)) + ['E(>%d)' % (row_max - 1)])
        fig_heatmap.update_xaxes(title_text="Allele 2", tickmode='array',
                                 tickvals=np.concatenate([np.array(range(0, column_max - 1, 5)), [column_max - 1]]),
                                 ticktext=list(range(0, column_max - 1, 5)) + ['E(>%d)' % (column_max - 1)])

        # generate motif
        generate_motif_report(output_dir, _key, motifs[_key], plots[_key], alignments[_key], fig_heatmap, fig_histogram, bgs)


if __name__ == '__main__':
    input_d, output_d, args = load_arguments()
    create_reports(input_d, output_d, args)
