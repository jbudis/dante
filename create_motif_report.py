import os
import argparse
import textwrap
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
  </tr>"""


def load_arguments():
    """
    Loads and parses arguments.
    :return: config_dict - config file in dictionary format.
             table_df - motif table in dataframe format.
             args - parsed arguments.
    """
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description=textwrap.dedent("""Python program to collect Dante report files and 
                                     create 'reversed' report for each motif, to make the comparison of multiple samples 
                                     easier"""))

    # add arguments
    parser.add_argument('input_dir', help='Path to directory with Dante reports')
    parser.add_argument('output_dir', help='Path to directory where the output will be stored')

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
                             motif_conf=c,reads_blue=reads_blue, reads_grey=reads_grey)


def generate_motif_report(path, key, samples):
    template = open('report/motif_report.html', 'r').read()
    rows = []
    key = key.replace('/', '-')

    for sample in samples:
        rows.append(generate_row(sample[0], sample[1], sample[2], sample[3], sample[4], sample[5], sample[6], sample[7]))

    with open('%s/report_%s.html' % (path, key), 'w') as f:
        table = motif_summary.format(table='\n'.join(rows))
        f.write(custom_format(template, motif=key.split('_')[0], seq=key.split('_')[1], motifs_content=table))


def create_reports(input_dir, output_dir):
    os.makedirs(output_dir, exist_ok=True)
    paths = []

    # find all all_profiles.txt files in root directory
    for root, dirs, files in os.walk(input_dir):
        for file in files:
            if file.endswith('.html'):
                paths.append(os.path.join(root, file))

    motifs = {}

    for path in paths:
        file = BeautifulSoup(open(path, 'r'), 'html.parser')
        fname = path.split('/')[-1].split('.')[0]

        for row in file.find(class_='tg').find_all('tr'):
            columns = row.find_all('td')

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

    for _key in motifs.keys():
        generate_motif_report(output_dir, _key, motifs[_key])


if __name__ == '__main__':
    input_d, output_d, args = load_arguments()
    create_reports(input_d, output_d)
