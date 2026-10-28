"""Render the research README to Word with editable equations and shared content.

Requires Python, python-docx and Pandoc (also distributed with Quarto).
Set PANDOC or pass --pandoc when Pandoc is not on PATH.
"""
from pathlib import Path
import argparse, hashlib, json, os, re, shutil, subprocess, tempfile
from docx import Document
from docx.shared import Inches, Pt, RGBColor
from docx.oxml import OxmlElement
from docx.oxml.ns import qn
from docx.enum.style import WD_STYLE_TYPE

ROOT=Path(__file__).resolve().parents[1]

def reference(path):
    d=Document();sec=d.sections[0]
    # Define styles referenced by Pandoc so Word and LibreOffice resolve them alike.
    for name in ['First Paragraph','Compact','Source Code']:
        if name not in d.styles:
            style=d.styles.add_style(name,WD_STYLE_TYPE.PARAGRAPH)
            style.base_style=d.styles['Normal']
    if 'Hyperlink' not in d.styles:
        style=d.styles.add_style('Hyperlink',WD_STYLE_TYPE.CHARACTER)
        style.font.color.rgb=RGBColor(0,70,130)
        style.font.underline=True
    sec.page_width=Inches(8.5);sec.page_height=Inches(11)
    sec.top_margin=sec.bottom_margin=Inches(.7)
    sec.left_margin=sec.right_margin=Inches(.72)
    for style in d.styles:
        if style.type==1:
            style.font.name='Calibri';style.font.size=Pt(11)
            style.font.color.rgb=RGBColor(0,0,0)
            style.paragraph_format.space_after=Pt(5)
            style.paragraph_format.line_spacing=1.04
    for name,size in [('Title',23),('Subtitle',12),('Heading 1',16),('Heading 2',13),('Heading 3',11)]:
        st=d.styles[name];st.font.size=Pt(size)
        st.paragraph_format.keep_with_next=True
        st.paragraph_format.space_before=Pt(10 if name.startswith('Heading') else 0)
    for e in list(d.styles.element.iter(qn('w:pBdr'))):e.getparent().remove(e)
    d.save(path)

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--pandoc',default=os.environ.get('PANDOC') or shutil.which('pandoc'))
    args=ap.parse_args()
    if not args.pandoc:ap.error('Pandoc is required: use --pandoc or set PANDOC.')
    source=(ROOT/'README.md').read_text(encoding='utf-8')
    title=source.splitlines()[0].removeprefix('# ')
    text='\n'.join(source.splitlines()[4:])
    text='\n'.join(line for line in text.splitlines() if not line.startswith('[Word manuscript]'))
    # No public repository exists yet; keep relative paths as explicit local references.
    text=re.sub(r'(?<!!)\[([^]]+)\]\((?!https?://|#)([^)]+)\)',r'\1 (\2)',text)
    # Standard Markdown plus TeX math lets Pandoc emit native Office Math.
    text=re.sub(r'!\[([^]]+)\]\(([^)]+)\)',r'![\1](\2){width=6.6in}',text)
    out=ROOT/'report/Bootstrap_LASSO_Manc_COJO_Report.docx'
    out.parent.mkdir(exist_ok=True)
    with tempfile.TemporaryDirectory(dir=ROOT/'work') as temp:
        temp=Path(temp);ref=temp/'reference.docx';reference(ref)
        md=temp/'manuscript.md';md.write_text(text,encoding='utf-8')
        subprocess.run([args.pandoc,str(md),'--from=markdown-implicit_figures+tex_math_dollars','--to=docx','--shift-heading-level-by=-1',
            '--resource-path='+str(ROOT),'--reference-doc='+str(ref),
            '--metadata=title:'+title,'--metadata=author:Vusumuzi Mabasa',
            '--metadata=subtitle:An experimental selection replacement evaluated with simulated phenotypes',
            '--output='+str(out)],cwd=ROOT,check=True)
    d=Document(out);d.core_properties.author='Vusumuzi Mabasa';d.core_properties.title=title
    for p in d.paragraphs:
        if p.style.name.startswith('Heading'):p.paragraph_format.keep_with_next=True
        if p._p.xpath('.//w:drawing'):p.paragraph_format.keep_with_next=True
        if p.text.startswith(('Figure ','Table ')):
            p.paragraph_format.space_after=Pt(8)
        if p.style.name=='Source Code':
            p.paragraph_format.line_spacing=1
            for run in p.runs:
                run.font.name='Consolas';run.font.size=Pt(9.5);run.font.color.rgb=RGBColor(0,0,0)
    # Rebuild Pandoc tables with complete native Word cell widths and standard styles.
    # This avoids importer-dependent placement of table text in LibreOffice.
    for old in list(d.tables):
        table=d.add_table(rows=len(old.rows), cols=len(old.columns))
        table.style='Table Grid'
        for src_row,dst_row in zip(old.rows,table.rows):
            for src,dst in zip(src_row.cells,dst_row.cells):
                dst.text=src.text
        old._tbl.addprevious(table._tbl)
        old._tbl.getparent().remove(old._tbl)
    for table in d.tables:
        table.autofit=True
        borders=OxmlElement('w:tblBorders')
        for side in ['top','left','bottom','right','insideH','insideV']:
            e=OxmlElement('w:'+side)
            for k,v in [('val','single'),('sz','4'),('color','D9D9D9')]:e.set(qn('w:'+k),v)
            borders.append(e)
        prior=table._tbl.tblPr.find(qn('w:tblBorders'))
        if prior is not None:table._tbl.tblPr.remove(prior)
        table._tbl.tblPr.insert_element_before(borders,'w:shd','w:tblLayout','w:tblCellMar','w:tblLook')
        for i,row in enumerate(table.rows):
            props=row._tr.get_or_add_trPr();props.append(OxmlElement('w:cantSplit'))
            if i==0 and props.find(qn('w:tblHeader')) is None:props.append(OxmlElement('w:tblHeader'))
            for cell in row.cells:
                cell.vertical_alignment=1
                if i==0:
                    shade=OxmlElement('w:shd');shade.set(qn('w:fill'),'E8EFF4');cell._tc.get_or_add_tcPr().insert_element_before(shade,'w:noWrap','w:tcMar','w:textDirection','w:tcFitText','w:vAlign')
                for p in cell.paragraphs:
                    if len(table.rows)<=10 or (len(table.rows)==17 and (i==0 or i%4!=0)):
                        p.paragraph_format.keep_with_next=(i==0)
                    p.paragraph_format.space_before=Pt(3);p.paragraph_format.space_after=Pt(3)
                    for r in p.runs:r.font.size=Pt(10);r.bold=(i==0)
    for e in list(d.styles.element.iter(qn('w:pBdr'))):e.getparent().remove(e)
    for section in d.sections:
        footer=section.footer.paragraphs[0];footer.alignment=2
        d.styles['Footer'].font.size=Pt(9)
        run=footer.add_run('Page ');run.font.size=Pt(9)
        field=OxmlElement('w:fldSimple');field.set(qn('w:instr'),'PAGE');footer._p.append(field)
    d.save(out)
    import docx
    (ROOT/'results/word_report_software.json').write_text(json.dumps({
        'builder':'src/build_word.py','source':'README.md',
        'source_sha256':hashlib.sha256(source.encode()).hexdigest(),
        'pandoc':subprocess.check_output([args.pandoc,'--version'],text=True).splitlines()[0],
        'python_docx':docx.__version__, 'author':'Vusumuzi Mabasa',
        'equations':'Native Office Math (OMML)',
        'analysis_validation':'results/audit.json'},indent=2)+'\n',encoding='utf-8')
    print(out)

if __name__=='__main__':main()
