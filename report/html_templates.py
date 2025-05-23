import gzip
import numpy as np
import re

contents = """
<table class="mtg" id="content-tg">
    <thead>
        <tr>
            <th class="mtg-s6z2">Motif</th>
        </tr>
    </thead>
    <tbody>
        {table}
    </tbody>
</table>
"""

make_datatable_string = """
<script>
    $(document).ready( function () {{
    $('#content-tg').DataTable();
}} );
</script>
"""

content_string = """ <tr>
    <td class="mtg-s6z2">
        <button class="tablinks" onclick="openTab(event, '{motif_name}'); $('.{motif_name}').trigger('content-change');">
            <a href="#{motif_name}">{motif}</a>
        </button>
    </td>
</tr>"""

content_string_empty = ""

row_string = """  <tr>
    <td class="tg-s6z2">{motif_name}</td>
    <td class="tg-s6z2">{str_seq}</td>
    <td class="tg-s6z2">{allele1}</td>
    <td class="tg-s6z2">{allele1_conf:.1%}</td>
    <td class="tg-s6z2">{allele2}</td>
    <td class="tg-s6z2">{allele2_conf:.1%}</td>
    <td class="tg-s6z2">{motif_conf:.1%}</td>
    <td class="tg-s6z2">{reads_blue}</td>
    <td class="tg-s6z2">{reads_grey}</td>
    <td class="tg-s6z2">{post_bases}</td>
    <td class="tg-s6z2">{post_reps}</td>
    <td class="tg-s6z2">{post_errors}</td>
    <td class="tg-s6z2">{motif_seq}</td>
  </tr>"""

row_string_empty_old = """  <tr>
    <td class="tg-s6z2">{motif_name}</td>
    <td class="tg-s6z2">{str_seq}</td>
    <td class="tg-s6z2">---</td>
    <td class="tg-s6z2">---</td>
    <td class="tg-s6z2">---</td>
    <td class="tg-s6z2">---</td>
    <td class="tg-s6z2">---</td>
    <td class="tg-s6z2">{reads_blue}</td>
    <td class="tg-s6z2">{reads_grey}</td>
    <td class="tg-s6z2">{post_bases}</td>
    <td class="tg-s6z2">{post_reps}</td>
    <td class="tg-s6z2">{post_errors}</td>
    <td class="tg-s6z2">{motif_seq}</td>
  </tr>"""

row_string_empty = ""

nomenclature_string = """
<tr>
    <td>{count}</td>
    <td>{ref}</td>
    {parts}
</tr>
"""

motif_summary = """
<div class="tabcontent" id="{motif_id}" style="display: none">
<h2 class="summary_nomenclatures">Nomenclatures</h2>
<table class="nomtg">
    <tbody>
        {nomenclatures}
    </tbody>
</table>
<h2 class="summary">Summary table</h2>
<table class="tg" id="tg-{motif_id}">
    <thead>
        <tr>
            <th class="tg-s6z2" rowspan="2">Motif</th>
            <th class="tg-s6z2" rowspan="2">STR<br>sequence</th>
            <th class="tg-s6z2" colspan="2">Allele 1</th>
            <th class="tg-s6z2" colspan="2">Allele 2</th>
            <th class="tg-s6z2" rowspan="2">Overall<br>confidence</th>
            <th class="tg-s6z2" colspan="2">Reads</th>
            <th class="tg-s6z2" colspan="3">Postfilter</th>
            <th class="tg-s6z2" rowspan="2">Sequence</th>
        </tr>
        <tr>
            <td class="tg-s6z2">prediction</td>
            <td class="tg-s6z2">confidence</td>
            <td class="tg-s6z2">prediction</td>
            <td class="tg-s6z2">confidence</td>
            <td class="tg-s6z2">full</td>
            <td class="tg-s6z2">partial</td>
            <td class="tg-s6z2">bases</td>
            <td class="tg-s6z2">modules</td>
            <td class="tg-s6z2">max. errors</td>
        </tr>
    </thead>
    <tbody>
        {table}
    </tbody>
</table>

<script>
    $(document).ready( function () {{
    $('#tg-{motif_id}').DataTable();
}} );
</script>
<p><a href="#content">Back to content</a></p>
{motifs}
</div>
"""

motif_summary_static = """
<div class="tabcontent" style="margin: 25px 0 25px 0">
<h2 class="summary">Summary table</h2>
<table class="tg">
    <thead>
        <tr>
            <th class="tg-s6z2" rowspan="2">Motif</th>
            <th class="tg-s6z2" rowspan="2">STR<br>sequence</th>
            <th class="tg-s6z2" colspan="2">Allele 1</th>
            <th class="tg-s6z2" colspan="2">Allele 2</th>
            <th class="tg-s6z2" rowspan="2">Overall<br>confidence</th>
            <th class="tg-s6z2" colspan="2">Reads</th>
            <th class="tg-s6z2" colspan="3">Postfilter</th>
            <th class="tg-s6z2" rowspan="2">Sequence</th>
        </tr>
        <tr>
            <td class="tg-s6z2">prediction</td>
            <td class="tg-s6z2">confidence</td>
            <td class="tg-s6z2">prediction</td>
            <td class="tg-s6z2">confidence</td>
            <td class="tg-s6z2">full</td>
            <td class="tg-s6z2">partial</td>
            <td class="tg-s6z2">bases</td>
            <td class="tg-s6z2">modules</td>
            <td class="tg-s6z2">max. errors</td>
        </tr>
    </thead>
    <tbody>
        {table}
    </tbody>
</table>
</div>
"""

motif_string = """<h2 id="{motif_name}">{motif}</h2>
<p>{sequence}</p>
postfilter: bases {post_bases} , repetitions {post_reps} , max. errors {errors}<br>
alleles: {result}<br>
<table class="plots">
    <tr>
        <td colspan="2">
            <img class="hist pic100" alt="{motif_name} repetitions" src="{motif_reps}" />
        </td>
        <td colspan="1">
            <img class="pcol pic100" alt="{motif_name} pcolor" src="{motif_pcolor}" />
        </td>
    </tr>
</table>
{alignment}
<p><a href="#content">Back to content</a></p>
"""

motif_stringb64 = """
<h2 id="{motif_name}">{motif}</h2>
<p>{sequence}</p>
postfilter: bases {post_bases} , repetitions {post_reps} , max. errors {errors}<br>
alleles: {result}<br>
<table class="plots">
    <tr>
        <td colspan="2">
            <div class="hist pic100 {motif_id}" id="hist-{motif_name}"></div>
            <script>
                {{
                    let hist_data = {motif_reps};
                    $(document).ready( function() {{
                        $('.{motif_id}').bind("content-change", function() {{
                            if (document.getElementById('{motif_id}').style.display === 'block') {{
                                Plotly.react('hist-{motif_name}', hist_data, {{}});
                            }}
                            else {{
                                Plotly.react('hist-{motif_name}', {{}}, {{}});
                            }}
                        }})
                    }})
                }}
            </script>
        </td>
        <td colspan="1">
            <div class="pcol pic100 {motif_id}" id="pcol-{motif_name}"></div>
            <script>
                {{
                    let pcol_data = {motif_pcolor};
                    $(document).ready( function() {{
                        $('.{motif_id}').bind("content-change", function() {{
                            if (document.getElementById('{motif_id}').style.display === 'block') {{
                                Plotly.react('pcol-{motif_name}', pcol_data, {{}});
                            }}
                            else {{
                                Plotly.react('pcol-{motif_name}', {{}}, {{}});
                            }}
                        }})
                    }})
                }}
            </script>
        </td>
    </tr>
</table>
{alignment}
<p><a href="#content">Back to content</a></p>
"""

motif_stringb64_reponly = """
<h2 id="{motif_name}">{motif}</h2>
<p>{sequence}</p>
postfilter: bases {post_bases} , repetitions {post_reps} , max. errors {errors}<br>
alleles: {result}<br>
<table class="plots">
    <tr>
        <td colspan="1">
            <div class="hist pic50 {motif_id}" id="hist2d-{motif_name}"></div>
            <script>
                {{
                    let hist2d_data = {motif_reps};
                    $(document).ready( function() {{
                        $('.{motif_id}').bind("content-change", function() {{
                            if (document.getElementById('{motif_id}').style.display === 'block') {{
                                Plotly.react('hist2d-{motif_name}', hist2d_data, {{}});
                            }}
                            else {{
                                Plotly.react('hist2d-{motif_name}', {{}}, {{}});
                            }}
                        }})
                    }})
                }}
            </script>
        </td>
    </tr>
</table>
{alignment}
<p><a href="#content">Back to content</a></p>
"""

motif_stringb64_static = """
<h2 id="{motif_name}">{motif}</h2>
<p>{sequence}</p>
postfilter: bases {post_bases} , repetitions {post_reps} , max. errors {errors}<br>
alleles: {result}<br>
<table class="plots">
    <tr>
        <td colspan="2">
            <div class="hist pic100 {motif_id}" id="hist-{motif_name}"></div>
            <script>
                Plotly.react('hist-{motif_name}', {motif_reps}, {{}});
            </script>
        </td>
        <td colspan="1">
            <div class="pcol pic100 {motif_id}" id="pcol-{motif_name}"></div>
            <script>
                Plotly.react('pcol-{motif_name}', {motif_pcolor}, {{}});
            </script>
        </td>
    </tr>
</table>
{alignment}
<p><a href="#content">Back to content</a></p>
"""

motif_stringb64_reponly_static = """
<h2 id="{motif_name}">{motif}</h2>
<p>{sequence}</p>
postfilter: bases {post_bases} , repetitions {post_reps} , max. errors {errors}<br>
alleles: {result}<br>
<table class="plots">
    <tr>
        <td colspan="1">
            <div class="hist pic50 {motif_id}" id="hist2d-{motif_name}"></div>
            <script>
                Plotly.react('hist2d-{motif_name}', {motif_reps}, {{}});
            </script>
        </td>
    </tr>
</table>
<p><a href="#content">Back to content</a></p>
"""

motif_string_empty_old = """<h2 id="{motif_name}">{motif}</h2>
<p>{sequence}</p>
NOT AVAILABLE
<p><a href="#content">Back to content</a></p>
"""

motif_string_empty = ""

alignment_string = """
  <p>{sequence}</p>
  {alignment}
  <hr>
"""

align_vis = """
  <details>
    <summary>{display_text}</summary>
    <div id="A{name}" class="align">press "Run with JS"</div>
    <script>
        var fasta = `{fasta}`;
        var seqs = msa.io.fasta.parse(fasta);
        var opts = {{
            el: document.getElementById("A{name}"),
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


def highlight_subpart(seq, highlight):
    """
    Highlights subpart of a motif sequence
    :param seq: str - motif sequence
    :param highlight: int/None/list(int) - part ot highlight
    :return: str, str - motif sequence with highlighted subpart, highlighted subpart
    """
    str_part = []
    if highlight is not None:
        highlight = np.array(highlight)
        split = seq.split(',')
        for h in highlight:
            str_part.append(split[h])
            split[h] = '<b><u>%s</u></b>' % split[h]
        return ','.join(split), ','.join(str_part)


def generate_row(motif, sequence, confidence, postfilter, reads_blue, reads_grey, highlight=None):
    """
    Generate rows of a summary table in html report.
    :param motif: str - motif name
    :param sequence: str - motif sequence
    :param confidence: tuple - motif confidences and allele predictions
    :param postfilter: dict - postfilter dict from config
    :param reads_blue: int - number of full reads
    :param reads_grey: int - number of partial reads
    :param highlight: list(int)/None - which part of seq to highlight
    :return: str - html string with rows of the summary table
    """
    sequence, subpart = highlight_subpart(sequence, highlight)

    # shorten sequence:
    keep = 10
    first = sequence.find(',')
    last = sequence.rfind(',')
    if first == -1:
        smaller_seq = sequence
    else:
        smaller_seq = '...' + sequence[first - keep:last + keep + 1] + '...'

    # errors:
    errors = postfilter['max_errors']
    if postfilter['max_errors_relative']:
        errors = '%.0f%%' % (postfilter['max_errors'] * 100)

    # fill templates:
    if confidence is None:
        return row_string_empty_old.format(post_bases=postfilter['bases'], post_reps=postfilter['repetitions'],
                                           motif_name=motif, motif_seq=smaller_seq, reads_blue=reads_blue,
                                           reads_grey=reads_grey, str_seq=subpart, post_errors=errors)
    else:
        (c, a1, a2, c1, c2, _, _, _, _) = confidence
        if a1 == 0 and a2 == 0:
            a1 = 'BG'
            a2 = 'BG'
        return row_string.format(post_bases=postfilter['bases'], post_reps=postfilter['repetitions'], motif_name=motif,
                                 motif_seq=smaller_seq, reads_blue=reads_blue, reads_grey=reads_grey, motif_conf=c,
                                 allele1=a1, allele2=a2, allele1_conf=c1, allele2_conf=c2, str_seq=subpart,
                                 post_errors=errors)


def get_alignment_name(alignment_file: str, allele: int) -> str:
    """
    Get alignment file of subpart of the alignment with allele count specified.
    :param alignment_file: str - alignment file name
    :param allele: int - allele repetition count
    :return: str - alignment file name for the subpart of alignment
    """
    # find where is .fasta
    fasta_index = alignment_file.rfind('.fasta')
    # insert '_aX' before .fasta
    return alignment_file[:fasta_index] + '_a' + str(allele) + alignment_file[fasta_index:]


def generate_motifb64(motif_name, description, sequence, repetition, pcolor, alignment, filtered_alignment, confidence,
                      postfilter, highlight=None, static=False, include_alignments=False):
    """
    Generate part of a html report for each motif.
    :param motif_name: str - motif name
    :param description: str - motif description
    :param sequence: str - motif sequence
    :param repetition: str - filename of repetitions figures
    :param pcolor: str - filename of pcolor figures
    :param alignment: str/None - filename of alignment file
    :param filtered_alignment: str/None - filename of filtered alignment file
    :param confidence: tuple - motif confidences and allele predictions
    :param postfilter: dict - postfilter dict from config
    :param highlight: list(int)/None - which part of seq to highlight
    :param static: bool - generate static HTML file?
    :param include_alignments: bool - add alignment logos to the main report?
    :return: (str, str) - content and main part of the html report for motifs
    """
    # prepare and generate alignments
    sequence, subpart = highlight_subpart(sequence, highlight)
    motif = '%s &ndash; %s' % (motif_name, description)
    motif_name = '%s_%s' % (motif_name.replace('/', '_'), ','.join(map(str, highlight)) if highlight is not None else 'mot')
    motif_clean = re.sub(r'[^\w_]', '', motif_name)
    align_html_a1 = ''
    align_html_a2 = ''
    if confidence is None:
        result = '-- (---.-%%) -- (---.-%%) total ---.-%%'
    else:
        (c, a1, a2, c1, c2, _, _, _, _) = confidence
        if (a1 == 'B' and a2 == 'B') or (a1 == 0 and a2 == 0):
            result = 'BG %5.1f%%' % (c * 100)
        else:
            result = '%2s (%5.1f%%) %2s (%5.1f%%) total %5.1f%%' % (str(a1), c1 * 100, str(a2), c2 * 100, c * 100)
            if alignment is not None:
                align_html_a1 = generate_alignment('%s_%s' % (motif_clean, str(a1)),
                                                   get_alignment_name(alignment, a1), motif_clean.split('_')[0],
                                                   'Allele 1 (%2s) alignment visualization' % str(a1))
                if a1 != a2:
                    align_html_a2 = generate_alignment('%s_%s' % (motif_clean, str(a2)),
                                                       get_alignment_name(alignment, a2), motif_clean.split('_')[0],
                                                       'Allele 2 (%2s) alignment visualization' % str(a2))

    # errors:
    errors = postfilter['max_errors']
    if postfilter['max_errors_relative']:
        errors = '%.0f%%' % (postfilter['max_errors'] * 100)

    # return content and picture parts:
    motif_templates = {'static': {'pcol': motif_stringb64_static, 'no-pcol': motif_stringb64_reponly_static},
                       'dynamic': {'pcol': motif_stringb64, 'no-pcol': motif_stringb64_reponly}}

    if repetition is not None:
        reps = open(repetition, 'r').read()
        align_html = generate_alignment(motif_clean, alignment, motif_clean.split('_')[0])
        filt_align_html = generate_alignment(motif_clean + '_filtered', filtered_alignment, motif_clean.split('_')[0],
                                             'Partial reads alignment visualization')
        # select template
        motif_template = motif_templates['static' if static else 'dynamic']['no-pcol' if pcolor is None else 'pcol']

        # read pcolor if available
        pcol = '' if pcolor is None else open(pcolor, 'r').read()

        # return filled valid template
        if include_alignments:
            return content_string.format(motif_name=motif_clean.rsplit('_', 1)[0], motif=motif), \
                motif_template.format(post_bases=postfilter['bases'], post_reps=postfilter['repetitions'],
                                      motif_name=motif_clean, motif_id=motif_clean.rsplit('_', 1)[0], motif=motif,
                                      motif_reps=reps, result=result, motif_pcolor=pcol,
                                      alignment=align_html + align_html_a1 + align_html_a2 + filt_align_html,
                                      sequence=sequence, errors=errors), \
                (motif, alignment_string.format(sequence=sequence,
                                                alignment=align_html + align_html_a1 + align_html_a2 + filt_align_html))
        else:
            return content_string.format(motif_name=motif_clean.rsplit('_', 1)[0], motif=motif), \
                motif_template.format(post_bases=postfilter['bases'], post_reps=postfilter['repetitions'],
                                      motif_name=motif_clean, motif_id=motif_clean.rsplit('_', 1)[0], motif=motif,
                                      motif_reps=reps, result=result, motif_pcolor=pcol,
                                      alignment=f'<p><a href="{motif_name}/alignments.html">Link to alignments</a></p>',
                                      sequence=sequence, errors=errors), \
                (motif, alignment_string.format(sequence=sequence,
                                                alignment=align_html + align_html_a1 + align_html_a2 + filt_align_html))
    else:
        return content_string_empty.format(motif_name=motif_clean.rsplit('_', 1)[0], motif=motif), \
            motif_string_empty.format(post_bases=postfilter['bases'], post_reps=postfilter['repetitions'],
                                      motif_name=motif_clean, motif=motif, sequence=sequence, errors=errors), \
            (motif, '')


def generate_alignment(motif: str, alignment_file: str, motif_id: str, display_text: str = 'Click to toggle alignment visualization'):
    """
    Generate HTML code for the fancy alignment.
    :param motif: str - name of the motif
    :param alignment_file: str - filename of the alignment file
    :param motif_id: str - motif identification
    :param display_text: str - string to display when the alignment is hidden
    :return: str - code of the fancy alignment
    """
    if alignment_file is None:
        return ''

    try:
        with gzip.open(alignment_file, 'rt') if alignment_file.endswith('.gz') else open(alignment_file) as f:
            string = f.read()
        debug = string.find('#')

        return align_vis.format(fasta=string[:debug], name=motif, motif_id=motif_id, display_text=display_text)
    except (IOError, TypeError, AttributeError):
        try:
            with gzip.open(alignment_file + '.gz', 'rt') as fgz:
                string = fgz.read()
            debug = string.find('#')

            return align_vis.format(fasta=string[:debug], name=motif, motif_id=motif_id, display_text=display_text)
        except (IOError, TypeError, AttributeError):
            return ''
