import base64
import numpy as np
import json
import os

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

motif_summary = """
<div class="tabcontent" id="{motif_name}">
<h2 id="summary">Summary table</h2>
<table class="tg" id="{motif_tg}_tg">
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
    $('#{motif_tg}_tg').DataTable();
}} );
</script>
<p><a href="#content">Back to content</a></p>
{motifs}
</div>
"""

motif_string = """<h2 id="{motif_name}">{motif}</h2>
{sequence}<br>
postfilter: bases {post_bases} , repetitions {post_reps} , max. errors {errors}<br>
alleles: {result}<br>
<img class="pic60" alt="{motif_name} repetitions" src="{motif_reps}" />
<img class="pic30" alt="{motif_name} pcolor" src="{motif_pcolor}" />
{alignment}
<p><a href="#content">Back to content</a></p>
"""

motif_stringb64 = """
<h2 id="{motif_name}">{motif}</h2>
{sequence}<br>
postfilter: bases {post_bases} , repetitions {post_reps} , max. errors {errors}<br>
alleles: {result}<br>
<div class="row">
<div class="pic60" id="plotly_{motif_name}"></div>
<script>
    Plotly.newPlot('plotly_{motif_name}', {motif_reps}, {{}});
</script>
<img class="pic30" alt="{motif_name} pcolor" src="data:image/png;base64,{motif_pcolor}" />
</div>
{alignment}
<p><a href="#content">Back to content</a></p>
"""

motif_stringb64_reponly = """
<h2 id="{motif_name}">{motif}</h2>
{sequence}<br>
postfilter: bases {post_bases} , repetitions {post_reps} , max. errors {errors}<br>
alleles: {result}<br>
<img class="pic30" alt="{motif_name} repetitions" src="data:image/png;base64,{motif_reps}" />
{alignment}
<p><a href="#content">Back to content</a></p>
"""

motif_string_empty_old = """<h2 id="{motif_name}">{motif}</h2>
{sequence}<br>
NOT AVAILABLE
<p><a href="#content">Back to content</a></p>
"""

motif_string_empty = ""

align_vis = """
  <details>
    <summary>{display_text}</summary>
    <div id="A{name}">press "Run with JS"</div>
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

content_string = """
<button class="tablinks" onclick="openTab(event, '{motif_name}')">
<a href="#{motif_name}">{motif}</a>
</button>"""

content_string_empty = ""


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


def generate_motifb64(motif_name, description, sequence, repetition, pcolor, alignment, confidence, postfilter, highlight=None):
    """
    Generate part of a html report for each motif.
    :param motif_name: str - motif name
    :param description: str - motif description
    :param sequence: str - motif sequence
    :param repetition: str - filename of repetitions figures
    :param pcolor: str - filename of pcolor figures
    :param alignment: str - filename of alignment file
    :param confidence: tuple - motif confidences and allele predictions
    :param postfilter: dict - postfilter dict from config
    :param highlight: list(int)/None - which part of seq to highlight
    :return: (str, str) - content and main part of the html report for motifs
    """
    # prepare and generate alignments
    sequence, subpart = highlight_subpart(sequence, highlight)
    motif = '%s &ndash; %s' % (motif_name, description)
    motif_name = '%s_%s' % (motif_name, ','.join(map(str, highlight)) if highlight is not None else 'mot')
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
            align_prefix = alignment[:alignment.rfind('.')]
            align_html_a1 = generate_alignment('%s_%s' % (motif_name, str(a1)), '%s_a%s.fasta' % (align_prefix, str(a1)), "Allele 1 (%2s) alignment visualization" % str(a1))
            if a1 != a2:
                align_html_a2 = generate_alignment('%s_%s' % (motif_name, str(a2)), '%s_a%s.fasta' % (align_prefix, str(a2)), "Allele 2 (%2s) alignment visualization" % str(a2))

    # errors:
    errors = postfilter['max_errors']
    if postfilter['max_errors_relative']:
        errors = '%.0f%%' % (postfilter['max_errors'] * 100)

    # return content and picture parts:
    if repetition is not None:
        if postfilter['index_rep2'] != 'no':
            reps = base64.b64encode(open(repetition, "rb").read())
            reps = reps.decode("utf-8")
        else:
            reps = open(repetition, 'r').read()
        align_html = generate_alignment(motif_name, alignment)

        if pcolor is not None:
            pcol = base64.b64encode(open(pcolor, "rb").read())
            pcol = pcol.decode("utf-8")
            return content_string.format(motif_name=motif_name.split('_')[0], motif=motif, sequence=sequence), \
                   motif_stringb64.format(post_bases=postfilter['bases'], post_reps=postfilter['repetitions'],
                                          motif_name=motif_name, motif=motif, motif_reps=reps, result=result,
                                          motif_pcolor=pcol, alignment=align_html + align_html_a1 + align_html_a2,
                                          sequence=sequence, errors=errors)
        else:
            return content_string.format(motif_name=motif_name.split('_')[0], motif=motif, sequence=sequence), \
                   motif_stringb64_reponly.format(post_bases=postfilter['bases'], post_reps=postfilter['repetitions'],
                                                  motif_name=motif_name, motif=motif, motif_reps=reps, result=result,
                                                  alignment=align_html + align_html_a1 + align_html_a2,
                                                  sequence=sequence, errors=errors)

    else:
        return content_string_empty.format(motif_name=motif_name, motif=motif, sequence=sequence), \
               motif_string_empty.format(post_bases=postfilter['bases'], post_reps=postfilter['repetitions'],
                                         motif_name=motif_name, motif=motif, sequence=sequence, errors=errors)


def generate_alignment(motif, alignment_file, display_text="Click to toggle alignment visualization"):
    """
    Generate HTML code for the fancy alignment.
    :param motif: str - name of the motif
    :param alignment_file: str - filename of the alignment file
    :param display_text: str - string to display when the alignment is hidden
    :return: str - code of the fancy alignment
    """
    try:
        with open(alignment_file) as f:
            string = f.read()
        debug = string.find('#')
        return align_vis.format(fasta=string[:debug], name=motif, display_text=display_text)
    except (IOError, TypeError):
        return ""
