"""Writes paper_protocol/Image_analysis_protocol.docx (methods, figure legends, references).

Numbers come from figures/stats_summary.csv, written by make_figures.py.
Run with the `axon` environment:  ~/axon/bin/python make_protocol_docx.py
"""
import os
import pandas as pd
from docx import Document
from docx.shared import Pt, Inches, RGBColor
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.oxml.ns import qn

OUT = '/Users/ctuna/Desktop/axon_codes/paper_protocol'
FIG = os.path.join(OUT, 'figures')
st = pd.read_csv(os.path.join(FIG, 'stats_summary.csv'))


def row(section, metric_start, comp_start):
    m = st[(st.section == section) & st.metric.str.startswith(metric_start) & st.comparison.str.startswith(comp_start)]
    assert len(m) == 1, (section, metric_start, comp_start, len(m))
    return m.iloc[0]


def fp(p):
    if p < 0.001:
        return 'P < 0.001'
    return f'P = {p:.3f}'.rstrip('0').rstrip('.') if p >= 0.01 else f'P = {p:.3f}'


br = {e: row('Axon breaks', 'confirmed breaks per 100', f'Exp. {e} ') for e in (1, 2, 3, 4)}
br_all = row('Axon breaks', 'confirmed breaks per 100 µm, BMI1/SCR', 'BMI1 vs SCR across')
pu = {e: row('Puncta', 'puncta per 1,000 µm²', f'Exp. {e} ') for e in (1, 2, 3, 4)}
pu_all = row('Puncta', 'puncta per 1,000 µm², BMI1/SCR', 'BMI1 vs SCR across')
mo = {k: row('Morphology (2WPT CZI)', k, 'BMI1 vs SCR same mouse')
      for k in ['dense_frac', 'iso_branch_points_per_100um', 'iso_diam_median_um', 'iso_pct_branches_swollen']}
cells = {k: st[(st.metric == 'cell bodies') & (st.comparison == k)].iloc[0]
         for k in ['ID06902 scr', 'ID06902 bmi1', 'ID06914 scr']}

doc = Document()
sec = doc.sections[0]
sec.page_width, sec.page_height = Inches(8.27), Inches(11.69)           # A4
for side in ('left_margin', 'right_margin', 'top_margin', 'bottom_margin'):
    setattr(sec, side, Inches(1))

base = doc.styles['Normal']
base.font.name = 'Arial'; base.font.size = Pt(10.5)
base.element.rPr.rFonts.set(qn('w:eastAsia'), 'Arial')
base.paragraph_format.space_after = Pt(6)
base.paragraph_format.line_spacing = 1.15
for lvl, size in ((1, 14), (2, 12), (3, 11)):
    hs = doc.styles[f'Heading {lvl}']
    hs.font.name = 'Arial'; hs.font.size = Pt(size); hs.font.bold = True
    hs.font.color.rgb = RGBColor(0, 0, 0)
    rf = hs.element.rPr.rFonts
    for att in ('w:asciiTheme', 'w:hAnsiTheme', 'w:eastAsiaTheme', 'w:cstheme'):
        rf.attrib.pop(qn(att), None)
    rf.set(qn('w:eastAsia'), 'Arial')
    hs.paragraph_format.space_before = Pt(14 if lvl == 1 else 10)
    hs.paragraph_format.space_after = Pt(4)
    hs.paragraph_format.keep_with_next = True


def para(text, italic=False, size=None, align=None, bold_lead=None):
    p = doc.add_paragraph()
    if bold_lead:
        r = p.add_run(bold_lead); r.bold = True
        if size: r.font.size = Pt(size)
    r = p.add_run(text); r.italic = italic
    if size: r.font.size = Pt(size)
    if align: p.alignment = align
    return p


def note(text):
    p = doc.add_paragraph()
    r = p.add_run(text); r.italic = True; r.font.color.rgb = RGBColor(0xA0, 0x2B, 0x20)
    return p


def bullets(items):
    for it in items:
        p = doc.add_paragraph(style='List Bullet')
        if isinstance(it, tuple):
            r = p.add_run(it[0]); r.bold = True
            p.add_run(it[1])
        else:
            p.add_run(it)


def figure(path, number, title, legend):
    doc.add_page_break()
    p = doc.add_paragraph(); p.alignment = WD_ALIGN_PARAGRAPH.CENTER
    p.add_run().add_picture(path, width=Inches(6.27))
    p = doc.add_paragraph()
    r = p.add_run(f'Figure {number}. {title} '); r.bold = True; r.font.size = Pt(9.5)
    r = p.add_run(legend); r.font.size = Pt(9.5)


# ── title ───────────────────────────────────────────────────────────────────
t = doc.add_paragraph(); t.alignment = WD_ALIGN_PARAGRAPH.LEFT
r = t.add_run('Image analysis: axon integrity, tubulin puncta and morphology of grafted neurons')
r.bold = True; r.font.size = Pt(16)
para('Methods text, figures and references for the manuscript. Draft for the authors.', italic=True)
note('Text in red italics marks details that the wet-lab authors should complete or confirm '
     '(antibodies, culture and grafting protocols, number of independent cultures). Remove these notes before submission.')

# ── 1 acquisition ───────────────────────────────────────────────────────────
doc.add_heading('1. Image acquisition and data organisation', 1)
para('All images were acquired on a Zeiss LSM 800 confocal microscope (Axio Observer Z1/7 stand) with a '
     'C Plan-Apochromat 63x/1.40 oil objective and saved as CZI files. Z-stacks were taken with a 1 µm step. '
     'Images were read in Python with pylibCZIrw. For every image, the channel of interest was combined into one '
     '2D image by maximum intensity projection over all Z planes. Channels were selected by their name in the '
     'file metadata, because the channel order changed between acquisitions.')
para('In vitro cultures. Neurons were stained for βIII-tubulin (AF647 channel) and DAPI. Three conditions were '
     'compared: BMI1 knockdown (BMI1), scrambled control (SCR) and non-infected neurons (NI). The analysis was '
     'repeated on four independent experiments (Exp. 1 to Exp. 4). Exp. 4 contains BMI1 and SCR only. '
     'The pixel size was 87.7 nm (Exp. 1), 70.6 nm (Exp. 3 and 4) and 35.3 to 87.7 nm (Exp. 2). '
     'Airyscan-processed images of Exp. 2 were analysed separately from the regular confocal images.')
note('Please add: antibodies and dilutions, culture and infection protocol, day of fixation, and how many '
     'coverslips or wells each experiment contains.')
para('Grafted neurons. Neurons grafted into the mouse brain were imaged two weeks after transplantation '
     '(AF594 label, DAPI). Mouse ID06902 received both BMI1 and SCR neurons; mouse ID06914 received SCR neurons only '
     'and serves as a second SCR reference. Regular confocal images (72 fields: 33 BMI1, 21 SCR in mouse ID06902; '
     '18 SCR in mouse ID06914) are shown here; Airyscan-processed single planes were analysed separately and gave '
     'the same conclusions.')
note('Please add: grafting procedure, tissue processing, section thickness and the AF594 label.')

# ── 2 axon breaks ───────────────────────────────────────────────────────────
doc.add_heading('2. Axon break analysis', 1)
para('We propose an automatic method that counts breaks along βIII-tubulin positive axons. It has three steps: '
     'a pretrained deep learning model finds the axons, the axon mask is reduced to centre lines, and pairs of '
     'axon ends that face each other across a dark gap are counted as breaks (Figure 1A).')
doc.add_heading('2.1 Axon segmentation', 2)
para('Each projected image was resampled to a pixel size of 87.7 nm, so that all distances below have the same '
     'physical length in every experiment. Axons were segmented with FAST, a pretrained U-Net++ network with a '
     'ResNet34 encoder that segments filamentous structures in fluorescence images (Aljapur et al., 2026). The '
     'model weights were downloaded from the BioImage Model Zoo (model "reliable-waterbuffalo"; Ouyang et al., '
     '2022). The model normalises the image itself, so no background subtraction was applied. The network returns '
     'five classes; the probabilities of the two filament classes (fibres and filopodia) were added, and pixels '
     'above 0.5 formed the axon mask.')
para('We chose FAST after comparing it on Exp. 1 with classical ridge filters used for neurites (Frangi et al., '
     '1998; Sato et al., 1998; Meijering et al., 2004) and with two other pretrained models (DeepD3, Fernholz et '
     'al., 2024; vesselFM, Wittmann et al., 2025). A Frangi filter with an Otsu threshold traced only about one '
     'sixth of the axon network that the other methods found, and most of its gaps fell on dim but continuous '
     'axons. FAST traced almost the whole network without tuning and gave the largest share of real gaps.')
doc.add_heading('2.2 Break detection', 2)
para('The axon mask was thinned to one-pixel centre lines (Zhang and Suen, 1984; scikit-image, van der Walt et al., '
     '2014). Pieces shorter than 1.3 µm were removed. Axon ends were the centre-line pixels with exactly one '
     'neighbour. Two ends from different pieces were paired as a break when all of the following held:')
bullets(['they were 0.9 to 7 µm apart;',
         'the local direction of each axon end, measured over the last 1.3 µm, pointed towards the other end within 40°;',
         'the straight line between them crossed almost no axon signal (at most 10% of its length), so a crossing axon was not taken as a break;',
         'each axon end was used in one break only, with the closest pairs taken first.'])
doc.add_heading('2.3 Confirmation on the raw image', 2)
para('A gap in the mask can also come from an axon that is only dim. Each break was therefore checked on the '
     'original image. We measured the intensity along the middle part of the gap, the intensity of the axon at '
     'both ends, and the background outside the axons. A break was kept (confirmed) when the gap was at most '
     'about one third as bright as the axon next to it (gap contrast below 0.35). The readout per image is the '
     'number of confirmed breaks per 100 µm of axon centre line.')
doc.add_heading('2.4 Background subtraction', 2)
para('No background subtraction was applied before segmentation or spot detection. Both FAST and Spotiflow '
     'normalise the image intensities themselves, and the break confirmation step already subtracts the local '
     'background. We checked this on Exp. 3 (54 images) by repeating both analyses after a white top-hat '
     'background subtraction (5 µm). Per-image values were almost unchanged (rank correlation with the original '
     'values 0.97 for confirmed breaks and 1.00 for puncta counts), and the BMI1 versus SCR result stayed the same '
     '(breaks 2.8 vs 3.0 fold, puncta 2.6 vs 2.6 fold, P < 0.01 in both cases).')

# ── 3 puncta ────────────────────────────────────────────────────────────────
doc.add_heading('3. Detection of tubulin puncta', 1)
para('Bright βIII-tubulin puncta along the axons were detected with Spotiflow, a deep learning spot detector, '
     'using its pretrained "general" model with default settings (Dominguez Mantes et al., 2025). The model was '
     'applied to the raw maximum intensity projection of the tubulin channel, without background subtraction. '
     'Each detected spot was counted once. Because the field size differs between experiments, counts are given '
     'per 1,000 µm² of imaged area (Figure 2A).')

# ── 4 morphology ────────────────────────────────────────────────────────────
doc.add_heading('4. Morphology of grafted neurons', 1)
para('We measured the shape of the AF594 labelled grafted neurons with an automatic pipeline. It finds the '
     'labelled cells and their extensions, then measures extension diameter, branching and swellings. The '
     'criteria follow the dendrite and axon measures used in an earlier study of BMI1 deficiency in human '
     'neurons (Flamier et al., 2018).')
doc.add_heading('4.1 Segmentation', 2)
bullets([('Resampling and background. ', 'Images were resampled to 70.6 nm pixels. A smooth background '
          '(Gaussian blur, 5 µm) was subtracted and the result was divided by the image noise, giving a '
          'signal-to-noise image.'),
         ('Labelled structures. ', 'The FAST filament mask (section 2.1) was combined with all pixels above '
          'three times the noise. A connected piece was kept only if it contained clearly bright signal (above six '
          'times the noise).'),
         ('Cell bodies. ', 'Cell bodies were found with Cellpose-SAM (Pachitariu et al., 2025; Stringer et al., '
          '2021), run at 0.12 µm per pixel. Objects of 30 to 900 µm² that were compact (solidity at least 0.85) and '
          'AF594 positive were kept.'),
         ('Extensions. ', 'The labelled mask without the cell bodies was thinned to centre lines and split into '
          'branches with skan (Nunez-Iglesias et al., 2018). End branches shorter than 2 µm were removed.')])
doc.add_heading('4.2 Measurements', 2)
bullets([('Dense clusters. ', 'Many BMI1 fields contained dense graft cores, where cells and processes are packed '
          'together. For every point of an extension, the fraction of labelled area in a 20 µm window was measured; '
          'above 0.2 the point was assigned to a dense cluster. The fraction of extension length inside dense '
          'clusters was reported per image. All shape measures below use only isolated extensions, outside dense '
          'clusters, because crossings inside a cluster would be counted as branch points.'),
         ('Branching. ', 'Number of branch points per 100 µm of extension.'),
         ('Extension diameter. ', 'Every 0.5 µm along each branch, the intensity profile across the branch was '
          'measured over 6 µm, and the diameter was taken as its full width at half maximum. This measure does not '
          'depend on a threshold, so a difference in brightness does not change it. The median per image is reported.'),
         ('Swellings. ', 'As in Flamier et al. (2018), we counted branches with at least one swelling. A swelling '
          'was a place where the branch was more than twice as wide as its own median diameter (and at least '
          '0.3 µm wider), brighter than 1.2 times its median intensity, wider than its neighbours, and at least '
          '2 µm away from a branch end.')])
para('Dendrite length and area per neuron, as traced in Flamier et al. (2018), require neurons whose whole arbor is '
     'visible. After visual inspection only two BMI1 neurons met this condition, so these per-neuron measures are '
     'not reported.')

# ── 5 statistics ────────────────────────────────────────────────────────────
doc.add_heading('5. Statistics', 1)
para('Each image (field) gave one value. Within an experiment, BMI1 and SCR were compared with a two-sided '
     'Mann-Whitney U test. To combine the four in vitro experiments, the BMI1 mean of each experiment was divided '
     'by its SCR mean, and the log2 ratios were tested against 0 with a one-sample t test (N = 4 experiments). '
     'For grafted neurons, BMI1 was compared with SCR from the same mouse and with SCR from the second mouse '
     '(Mann-Whitney U test). Plots show each image as a dot with the mean ± SEM. Significance: * P < 0.05, '
     '** P < 0.01, *** P < 0.001, ns not significant. Analyses used Python 3.13, NumPy, pandas and SciPy '
     '(Virtanen et al., 2020).')

# ── 6 results ───────────────────────────────────────────────────────────────
doc.add_heading('6. Results', 1)
doc.add_heading('BMI1 knockdown increases axon breaks', 2)
para(f'BMI1 neurons had more confirmed axon breaks than SCR neurons in all four experiments (Figure 1B). The '
     f'difference was significant in Exp. 2 ({fp(br[2].p)}), Exp. 3 ({fp(br[3].p)}) and Exp. 4 ({fp(br[4].p)}), '
     f'and showed the same direction in Exp. 1, which had only four SCR fields ({fp(br[1].p)}). Across the four '
     f'experiments, BMI1 neurons had {br_all.ratio:.1f} times more breaks than SCR neurons ({fp(br_all.p)}; '
     f'Figure 1C). Non-infected neurons were close to SCR. Airyscan images of Exp. 2 gave the same result.')
doc.add_heading('BMI1 knockdown increases tubulin puncta', 2)
para(f'BMI1 neurons showed more βIII-tubulin puncta than SCR neurons in all four experiments (Figure 2B; Exp. 1 '
     f'{fp(pu[1].p)}, Exp. 2 {fp(pu[2].p)}, Exp. 3 {fp(pu[3].p)}, Exp. 4 {fp(pu[4].p)}). The increase was very '
     f'consistent, {pu_all.ratio:.1f} fold on average ({fp(pu_all.p)}; Figure 2C). In the Airyscan-processed '
     f'images of Exp. 2, puncta counts did not differ between conditions.')
doc.add_heading('Grafted BMI1 neurons form dense clusters and branch more', 2)
para(f'In the brain, BMI1 neurons were mostly found in dense clusters, while SCR neurons were isolated (Figure 3A). '
     f'A median of {mo["dense_frac"].median_bmi * 100:.0f}% of the BMI1 extension length lay in dense clusters, '
     f'against {mo["dense_frac"].median_scr * 100:.0f}% and {mo["dense_frac"].median_scr2 * 100:.0f}% for SCR in '
     f'the two mice ({fp(mo["dense_frac"].p)} and {fp(mo["dense_frac"].p_vs_scr2)}; Figure 3B). '
     f'{int(cells["ID06902 bmi1"].n_in_cluster)} of {int(cells["ID06902 bmi1"].n_somata)} BMI1 cell bodies sat '
     f'inside a cluster. Isolated BMI1 extensions had more branch points than SCR extensions from the same mouse '
     f'({mo["iso_branch_points_per_100um"].median_bmi:.1f} vs {mo["iso_branch_points_per_100um"].median_scr:.1f} '
     f'per 100 µm, {fp(mo["iso_branch_points_per_100um"].p)}; Figure 3C). BMI1 extensions were also thicker than '
     f'SCR extensions in the same mouse ({mo["iso_diam_median_um"].median_bmi:.2f} vs '
     f'{mo["iso_diam_median_um"].median_scr:.2f} µm, {fp(mo["iso_diam_median_um"].p)}), but SCR extensions in the '
     f'second mouse were as thick ({mo["iso_diam_median_um"].median_scr2:.2f} µm; Figure 3D). Swellings were rare '
     f'and did not differ between groups (Figure 3E). All BMI1 neurons come from one mouse, so these in vivo results '
     f'should be confirmed in more animals.')

# ── figures ─────────────────────────────────────────────────────────────────
figure(os.path.join(FIG, 'Figure1_axon_breaks.png'), 1, 'BMI1 knockdown increases axon breaks in cultured neurons.',
       '(A) βIII-tubulin staining of SCR and BMI1 neurons (Exp. 3; maximum intensity projections of fields close '
       'to the group median). Left, whole field; yellow circles mark confirmed breaks. Middle, the dashed square '
       'enlarged. Right, the same area with the detected axon centre lines: each colour is one continuous axon '
       'segment, and white circles mark confirmed breaks. SCR axons form long continuous segments, while BMI1 '
       'axons are beaded and broken into many short pieces. Scale bars, 20 µm (field) and 5 µm (zoom). (B) Confirmed breaks per 100 µm of axon in each experiment. Each dot is one image; lines '
       f'show the mean ± SEM. Exp. 1: n = {int(br[1].n_scr)}, {int(br[1].n_ni)} and {int(br[1].n_bmi)} images '
       f'(SCR, NI, BMI1); Exp. 2: n = {int(br[2].n_scr)}, {int(br[2].n_ni)}, {int(br[2].n_bmi)}; Exp. 3: '
       f'n = {int(br[3].n_scr)}, {int(br[3].n_ni)}, {int(br[3].n_bmi)}; Exp. 4: n = {int(br[4].n_scr)} and '
       f'{int(br[4].n_bmi)} (SCR, BMI1). Two-sided Mann-Whitney U test, BMI1 vs SCR. (C) Mean of each experiment '
       'relative to its SCR mean (symbols, one per experiment); bars show the mean ± SEM of the four experiments. '
       'One-sample t test on log2 ratios, N = 4. * P < 0.05, ** P < 0.01, *** P < 0.001, ns not significant.')
figure(os.path.join(FIG, 'Figure2_puncta.png'), 2, 'BMI1 knockdown increases βIII-tubulin puncta.',
       '(A) βIII-tubulin staining of SCR and BMI1 neurons (Exp. 3; fields close to the group median). Left, whole '
       'field. Middle, the dashed square enlarged. Right, the same area with the puncta detected by Spotiflow '
       '(yellow circles). SCR axons are smooth with few puncta; BMI1 axons carry rows of bright puncta. Scale bars, '
       '20 µm (field) and 2 µm (zoom). (B) Puncta per 1,000 µm² in each experiment. Each dot is one image; lines show the mean ± '
       'SEM. Same images as in Figure 1. Two-sided Mann-Whitney U test, BMI1 vs SCR. (C) Mean of each experiment '
       'relative to its SCR mean; bars show the mean ± SEM of the four experiments. One-sample t test on log2 '
       'ratios, N = 4. * P < 0.05, ** P < 0.01, *** P < 0.001.')
figure(os.path.join(FIG, 'Figure3_morphology.png'), 3, 'Morphology of grafted BMI1 and SCR neurons in the mouse brain.',
       '(A) AF594 labelled grafted neurons, two weeks after transplantation. Left, an isolated SCR neuron; middle, '
       'an isolated BMI1 neuron; right, a dense cluster of BMI1 neurons. Cyan outlines, cell bodies detected by '
       'Cellpose-SAM; yellow dashed outline, dense cluster. Scale bars, 10 µm (left, middle) and 50 µm (right). '
       '(B) Fraction of extension length inside dense clusters. (C to E) Measures of isolated extensions only: '
       'branch points per 100 µm (C), extension diameter measured as the full width at half maximum of the '
       'intensity profile (D), and percentage of branches with at least one swelling (E). Each dot is one image; '
       f'lines show the mean ± SEM. Mouse 1 (ID06902): n = {int(mo["dense_frac"].n_scr)} SCR and '
       f'{int(mo["dense_frac"].n_bmi)} BMI1 images; mouse 2 (ID06914): n = {int(mo["dense_frac"].n_scr2)} SCR images. '
       'Images without enough isolated extensions are not shown in C to E. Two-sided Mann-Whitney U test, BMI1 vs '
       'SCR of the same mouse and BMI1 vs SCR of the second mouse. * P < 0.05, *** P < 0.001, ns not significant.')

# ── software / references ───────────────────────────────────────────────────
doc.add_page_break()
doc.add_heading('Software and code availability', 1)
para('Python 3.13 with pylibCZIrw 6.0, NumPy 2.4, SciPy 1.17, scikit-image 0.26, pandas 3.0, PyTorch 2.12, '
     'Spotiflow 0.6.4, Cellpose 4.2 (Cellpose-SAM), skan 0.13 and matplotlib 3.10. FAST model weights: BioImage '
     'Model Zoo, "reliable-waterbuffalo" (CC BY 4.0). Analysis notebooks and the script that draws the figures '
     'are available from the authors.')
note('Add the repository link (for example GitHub or Zenodo) if the code is released with the paper.')

doc.add_heading('References', 1)
refs = [
    'Aljapur V, Gardner A, Carayanniotis J, Harris AR (2026). FAST: filamentous actin segmentation tool for quantifying cytoskeletal organization. J Cell Sci 139(4): jcs264265. https://doi.org/10.1242/jcs.264265',
    'Dominguez Mantes A, Herrera A, Khven I, Schlaeppi A, Kyriacou E, Tsissios G, et al. (2025). Spotiflow: accurate and efficient spot detection for fluorescence microscopy with deep stereographic flow regression. Nat Methods 22(7): 1495–1504. https://doi.org/10.1038/s41592-025-02662-x',
    'Fernholz MHP, Guggiana Nilo DA, Bonhoeffer T, Kist AM (2024). DeepD3, an open framework for automated quantification of dendritic spines. PLoS Comput Biol 20(2): e1011774. https://doi.org/10.1371/journal.pcbi.1011774',
    'Flamier A, El Hajjar J, Adjaye J, Fernandes KJ, Abdouh M, Bernier G (2018). Modeling late-onset sporadic Alzheimer’s disease through BMI1 deficiency. Cell Rep 23(9): 2653–2666. https://doi.org/10.1016/j.celrep.2018.04.097',
    'Frangi AF, Niessen WJ, Vincken KL, Viergever MA (1998). Multiscale vessel enhancement filtering. In: Medical Image Computing and Computer-Assisted Intervention (MICCAI 1998), Lecture Notes in Computer Science 1496: 130–137. https://doi.org/10.1007/BFb0056195',
    'Meijering E, Jacob M, Sarria JCF, Steiner P, Hirling H, Unser M (2004). Design and validation of a tool for neurite tracing and analysis in fluorescence microscopy images. Cytometry A 58(2): 167–176. https://doi.org/10.1002/cyto.a.20022',
    'Nunez-Iglesias J, Blanch AJ, Looker O, Dixon MW, Tilley L (2018). A new Python library to analyse skeleton images confirms malaria parasite remodelling of the red blood cell membrane skeleton. PeerJ 6: e4312. https://doi.org/10.7717/peerj.4312',
    'Ouyang W, Beuttenmueller F, Gómez-de-Mariscal E, Pape C, Burke T, Garcia-López-de-Haro C, et al. (2022). BioImage Model Zoo: a community-driven resource for accessible deep learning in bioimage analysis. bioRxiv. https://doi.org/10.1101/2022.06.07.495102',
    'Pachitariu M, Rariden M, Stringer C (2025). Cellpose-SAM: superhuman generalization for cellular segmentation. bioRxiv. https://doi.org/10.1101/2025.04.28.651001',
    'Sato Y, Nakajima S, Shiraga N, Atsumi H, Yoshida S, Koller T, et al. (1998). Three-dimensional multi-scale line filter for segmentation and visualization of curvilinear structures in medical images. Med Image Anal 2(2): 143–168. https://doi.org/10.1016/S1361-8415(98)80009-1',
    'Stringer C, Wang T, Michaelos M, Pachitariu M (2021). Cellpose: a generalist algorithm for cellular segmentation. Nat Methods 18(1): 100–106. https://doi.org/10.1038/s41592-020-01018-x',
    'van der Walt S, Schönberger JL, Nunez-Iglesias J, Boulogne F, Warner JD, Yager N, et al. (2014). scikit-image: image processing in Python. PeerJ 2: e453. https://doi.org/10.7717/peerj.453',
    'Virtanen P, Gommers R, Oliphant TE, Haberland M, Reddy T, Cournapeau D, et al. (2020). SciPy 1.0: fundamental algorithms for scientific computing in Python. Nat Methods 17(3): 261–272. https://doi.org/10.1038/s41592-019-0686-2',
    'Wittmann B, Wattenberg Y, Amiranashvili T, Shit S, Menze B (2025). vesselFM: a foundation model for universal 3D blood vessel segmentation. Proceedings of the IEEE/CVF Conference on Computer Vision and Pattern Recognition (CVPR). arXiv: 2411.17386',
    'Zhang TY, Suen CY (1984). A fast parallel algorithm for thinning digital patterns. Commun ACM 27(3): 236–239. https://doi.org/10.1145/357994.358023',
]
for ref in refs:
    p = doc.add_paragraph(ref)
    p.paragraph_format.left_indent = Inches(0.3); p.paragraph_format.first_line_indent = Inches(-0.3)
    for r in p.runs:
        r.font.size = Pt(9.5)

out = os.path.join(OUT, 'Image_analysis_protocol.docx')
doc.save(out)
print('saved', out)
