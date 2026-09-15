INLINE=[]
NCOL=[0]
import re, os
from PIL import Image
doc='/home/ocs/matlab/Lizoo/doc/'
md=open(doc+'report_v13.md').read()

def md2html(t):
    out=[]; intable=False; para=[]; inli=False
    def flush():
        # join consecutive prose lines into ONE paragraph so the text reflows
        nonlocal inli
        if para:
            out.append('<p>'+' '.join(para)+'</p>'); para.clear()
        inli=False
    for line in t.split('\n'):
        raw=line.rstrip().replace('<','&lt;').replace('>','&gt;'); s=raw.strip()
        if s.startswith('@@FIG '):
            flush()
            body=s[6:]
            fn,_,cap=body.partition('|')
            fn=fn.strip(); cap=cap.strip()
            INLINE.append(fn)
            w,h=Image.open(doc+fn).size; hmm=W*h/w
            out.append(f'<div class="figpage"><img src="{doc+fn}" '
                       f'style="width:{W:.1f}mm;height:{hmm:.1f}mm">'
                       f'<p class="cap">{cap}<br><span class="fn">{fn}</span></p></div>')
            continue
        # a continuation of the bullet above: fold it back into that item
        if inli and raw.startswith(('  ','\t')) and s and not s.startswith(('- ','|','#')):
            out[-1]=out[-1][:-5]+' '+s+'</li>'
            continue
        if s.startswith('|'):
            flush()
            cells=[c.strip() for c in s.strip('|').split('|')]
            if set(''.join(cells))<=set('-: '): continue
            # A ragged table makes LibreOffice's layout spin forever, so every
            # row is padded or trimmed to the width of the header row.
            if not intable: NCOL[0]=len(cells)
            if len(cells)<NCOL[0]: cells=cells+['']*(NCOL[0]-len(cells))
            elif len(cells)>NCOL[0]: cells=cells[:NCOL[0]-1]+[' '.join(cells[NCOL[0]-1:])]
            tag='td'
            if not intable: out.append('<table width="100%" cellspacing="0">'); intable=True; tag='th'
            out.append('<tr>'+''.join(f'<{tag} nowrap>{c}</{tag}>' for c in cells)+'</tr>'); continue
        if intable: out.append('</table>'); intable=False
        if s.startswith('### '):  flush(); out.append(f'<h3>{s[4:]}</h3>')
        elif s.startswith('## '): flush(); out.append(f'<h2>{s[3:]}</h2>')
        elif s.startswith('# '):  flush(); out.append(f'<h1>{s[2:]}</h1>')
        elif s.startswith('---'): flush(); out.append('<hr>')
        elif s.startswith('- '):  flush(); out.append(f'<li>{s[2:]}</li>'); inli=True
        elif s=='':               flush()
        else:                     para.append(s)
    flush()
    if intable: out.append('</table>')
    h='\n'.join(out)
    h=re.sub(r'\*\*(.+?)\*\*', r'<b>\1</b>', h)
    h=re.sub(r'`(.+?)`', r'<code>\1</code>', h)
    h=re.sub(r'(?<!\*)\*([^*\n]+?)\*(?!\*)', r'<i>\1</i>', h)
    h=re.sub(r'((?:<li>.*?</li>\n?)+)', r'<ul>\1</ul>', h)
    return h

W=163.0
figs=[('report_v13_RMS.png',
   'Residual RMS against the OGLE catalogue magnitude, per axis and per field, after position and '
   'proper motion have been removed. Blue: all 287 / 253 calibration stars, every one freely fitted. '
   'Green: the subset with a Gaia RUWE below 1.4, used only for the tie. Red: the target. '
   'Black: running median.'),
  ('report_v13_pm_vs_gaia.png',
   'Our absolute proper motion against Gaia&rsquo;s, per axis and per field, after the affine tie with '
   'the <b>full six-parameter gauge</b> removed. Dashed line is 1:1. The scatter quoted in each panel '
   'is the robust standard deviation of the difference; it improves by a factor of three to four '
   'over a tie that removes only the constant. The target is the red star.'),
  ('report_v13_chi2.png',
   '&chi;<sup>2</sup> of the sidereal-month binned residuals for every calibration star, both axes '
   'combined, <b>not divided by DoF</b>. Each bin contributes (mean / standard error)&sup2; per axis '
   'and DoF = 2N<sub>bins</sub> &minus; 4. The red line marks the target, which sits in the body of '
   'the distribution rather than the tail.')]
RADEC=[('target','the target, I = 18.13'),
       ('m16_d24','comparison star I = 16.68, 24 pix from the target &mdash; the best measured of the set'),
       ('m16_d30','comparison star I = 16.63, 30 pix from the target')]
for F in ['BLG41','BLG01']:
    for tag,lab in RADEC:
        figs.append((f'report_v13_radec_{F}_{tag}.png',
            f'<b>{F}</b> &mdash; {lab}. <b>Position and proper motion removed</b>: these are the '
            f'residuals. Left: RA residual against time. Middle: Dec residual against time. Right: '
            f'the two against each other with time colour-coded, where a well-behaved star shows an '
            f'unordered cloud about the origin and any systematic drift would appear as an ordered '
            f'progression of colour; consecutive months are joined by a thin grey line. Sidereal-month bins; error bars are the standard error within '
            f'the bin scaled by the empirical factor &radic;(&chi;&sup2;/DoF) of the bin means, quoted in each panel title, so that they carry '
            f'the real bin-to-bin scatter (section 7). Bins with fewer than 10 epochs and epochs at '
            f'sec z &gt; 1.3 are dropped: the season-edge months, with a handful of high-airmass frames, otherwise dominate the plot.'))
ROLE=[('target','the target, I = 18.13'),
      ('m18_d04','I = 18.14, 4.5 pix away &mdash; <b>blended with the target</b>, see section 7'),
      ('m18_d17','I = 18.09, 17 pix away &mdash; poorly measured, see section 7'),
      ('m18_d48','I = 18.05, 48 pix away'),
      ('m16_d24','I = 16.68, 24 pix away'),
      ('m16_d30','I = 16.63, 30 pix away'),
      ('m16_d31','I = 16.75, 31 pix away')]
for F in ['BLG41','BLG01']:
    for tag,lab in ROLE:
        figs.append((f'report_v13_motion_{F}_{tag}.png',
            f'<b>{F}</b> &mdash; {lab}. Binned in one sidereal month, error bars are the standard error within '
            f'the bin scaled by the empirical factor &radic;(&chi;&sup2;/DoF) of the bin means (quoted in the panel titles, computed on the '
            f'proper-motion-removed residuals and applied to both rows); bins with fewer than 10 epochs and epochs at sec z &gt; 1.3 are dropped. Columns are the two axes, <b>X left and Y right</b>; the top row keeps the '
            f'proper motion, the bottom row removes it. The heading gives the absolute sky motion '
            f'with the full gauge removed; the legend gives it along the pixel axis, and <b>pixel X '
            f'runs opposite to RA</b>, so the two carry opposite signs in X.'))
_html=md2html(md)
figs=[x for x in figs if x[0] not in INLINE]
body=_html+'\n<hr><h2>Figures</h2>\n'
for i,(f,cap) in enumerate(figs,1):
    w,h=Image.open(doc+f).size; hmm=W*h/w
    body+=(f'<div class="fig"><img src="{doc+f}" style="width:{W:.1f}mm;height:{hmm:.1f}mm">'
           f'<p class="cap">Figure {i}. {cap}<br><span class="fn">{f}</span></p></div>\n')
css=f"""@page{{size:A4 portrait;margin:15mm}}
body{{font-family:Helvetica,Arial,sans-serif;font-size:9.5pt;line-height:1.30;text-align:justify}}
h1{{font-size:16pt;margin-bottom:2pt}}h2{{font-size:12pt;margin-top:13pt;border-bottom:1px solid #bbb}}
h3{{font-size:10.5pt;margin-top:10pt}}p{{margin:4pt 0}}
table{{border-collapse:collapse;margin:5pt 0;font-size:9pt;width:100%}}
td,th{{border:1px solid #999;padding:2pt 7pt;white-space:nowrap}}th{{background:#eee}}
code{{font-family:monospace;font-size:8.5pt;background:#f2f2f2;padding:0 2px}}
ul{{margin:3pt 0 3pt 0}}li{{margin-left:4pt}}
.fig{{page-break-inside:avoid;margin-top:6pt}}
.figpage{{page-break-inside:avoid;margin-top:0}}
.cap{{font-size:8.5pt;color:#333;margin-top:1pt;text-align:left}}
.fn{{font-family:monospace;font-size:7.5pt;color:#777}}"""
open(os.path.join(os.path.dirname(os.path.abspath(__file__)),'report_v13.html'),'w').write(
    f'<html><head><meta charset="utf-8"><style>{css}</style></head><body>{body}</body></html>')
print('v13 html rebuilt; paragraphs joined for reflow')
