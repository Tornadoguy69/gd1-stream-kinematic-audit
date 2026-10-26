import os
import sys
from reportlab.lib.pagesizes import letter
from reportlab.lib import colors
from reportlab.platypus import (
    SimpleDocTemplate, Paragraph, Spacer, Table, TableStyle, PageBreak, Image, KeepTogether
)
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.pdfgen import canvas

PDF_FILENAME = "D:/ARGON/gd1_stream_kinematic_audit.pdf"

class NumberedCanvas(canvas.Canvas):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._saved_page_states = []

    def showPage(self):
        self._saved_page_states.append(dict(self.__dict__))
        self._startPage()

    def save(self):
        num_pages = len(self._saved_page_states)
        for state in self._saved_page_states:
            self.__dict__.update(state)
            self.draw_page_decorations(num_pages)
            super().showPage()
        super().save()

    def draw_page_decorations(self, page_count):
        self.saveState()
        self.setFont("Helvetica", 8)
        self.setFillColor(colors.HexColor("#4A5568"))
        
        # Header (Pages 2+)
        if self._pageNumber > 1:
            self.drawString(54, letter[1] - 32, "The GD-1 Stellar Stream Kinematic Tail Audit")
            self.drawRightString(letter[0] - 54, letter[1] - 32, "Tornado guy (2026)")
            self.setStrokeColor(colors.HexColor("#CBD5E0"))
            self.setLineWidth(0.5)
            self.line(54, letter[1] - 36, letter[0] - 54, letter[1] - 36)

        # Footer
        footer_text = "Preprint draft, not peer reviewed"
        page_text = f"Page {self._pageNumber} of {page_count}"
        self.drawString(54, 32, footer_text)
        self.drawRightString(letter[0] - 54, 32, page_text)
        self.setStrokeColor(colors.HexColor("#CBD5E0"))
        self.setLineWidth(0.5)
        self.line(54, 40, letter[0] - 54, 40)
        self.restoreState()

def build_pdf():
    doc = SimpleDocTemplate(
        PDF_FILENAME,
        pagesize=letter,
        leftMargin=50,
        rightMargin=50,
        topMargin=42,
        bottomMargin=42
    )
    
    styles = getSampleStyleSheet()
    
    title_style = ParagraphStyle(
        'DocTitle',
        parent=styles['Normal'],
        fontName='Helvetica-Bold',
        fontSize=16,
        leading=20,
        textColor=colors.HexColor("#1A202C"),
        alignment=1, # Center
        spaceAfter=8
    )
    
    author_style = ParagraphStyle(
        'AuthorStyle',
        parent=styles['Normal'],
        fontName='Helvetica-Bold',
        fontSize=10,
        leading=14,
        textColor=colors.HexColor("#2B6CB0"),
        alignment=1,
        spaceAfter=12
    )
    
    abstract_heading = ParagraphStyle(
        'AbstractHeading',
        parent=styles['Normal'],
        fontName='Helvetica-Bold',
        fontSize=9.5,
        leading=12,
        textColor=colors.HexColor("#1A202C"),
        alignment=1,
        spaceAfter=4
    )

    abstract_style = ParagraphStyle(
        'AbstractBody',
        parent=styles['Normal'],
        fontName='Helvetica-Oblique',
        fontSize=7.4,
        leading=9.6,
        textColor=colors.HexColor("#2D3748"),
        leftIndent=16,
        rightIndent=16,
        spaceAfter=8
    )
    
    h1_style = ParagraphStyle(
        'Heading1_Custom',
        parent=styles['Normal'],
        fontName='Helvetica-Bold',
        fontSize=10.5,
        leading=13,
        textColor=colors.HexColor("#1A365D"),
        spaceBefore=7,
        spaceAfter=3,
        keepWithNext=True
    )

    h2_style = ParagraphStyle(
        'Heading2_Custom',
        parent=styles['Normal'],
        fontName='Helvetica-Bold',
        fontSize=9.0,
        leading=11.5,
        textColor=colors.HexColor("#2C5282"),
        spaceBefore=6,
        spaceAfter=2,
        keepWithNext=True
    )

    body_style = ParagraphStyle(
        'Body_Custom',
        parent=styles['Normal'],
        fontName='Helvetica',
        fontSize=7.8,
        leading=10.0,
        textColor=colors.HexColor("#2D3748"),
        spaceAfter=4
    )

    bullet_style = ParagraphStyle(
        'Bullet_Custom',
        parent=styles['Normal'],
        fontName='Helvetica',
        fontSize=7.6,
        leading=9.8,
        textColor=colors.HexColor("#2D3748"),
        leftIndent=12,
        spaceAfter=2
    )
    
    caption_style = ParagraphStyle(
        'CaptionStyle',
        parent=styles['Normal'],
        fontName='Helvetica',
        fontSize=7.2,
        leading=9.5,
        textColor=colors.HexColor("#4A5568"),
        spaceAfter=8
    )

    table_cell = ParagraphStyle(
        'TableCell',
        parent=styles['Normal'],
        fontName='Helvetica',
        fontSize=6.8,
        leading=8.5,
        textColor=colors.HexColor("#1A202C")
    )
    
    table_header = ParagraphStyle(
        'TableHeader',
        parent=styles['Normal'],
        fontName='Helvetica-Bold',
        fontSize=6.8,
        leading=8.5,
        textColor=colors.white
    )

    story = []
    
    # -------------------------------------------------------------------------
    # PAGE 1: Title, Abstract, Section 1, Section 2, Section 3
    # -------------------------------------------------------------------------
    story.append(Paragraph("The GD-1 Stellar Stream Kinematic Tail: Threshold Fragility, Spectroscopic Binarity, and Kinematics across the Spur in DESI DR2", title_style))
    story.append(Paragraph("Tornado guy", author_style))
    
    affiliation_style = ParagraphStyle(
        'AffiliationStyle',
        parent=styles['Normal'],
        fontName='Helvetica-Oblique',
        fontSize=8.0,
        leading=11.0,
        textColor=colors.HexColor("#4A5568"),
        alignment=1,
        spaceAfter=10
    )
    story.append(Paragraph("Independent Astrophysical Research<br/><i>Kinematic audit, Bayesian profile modeling &amp; pipeline orchestration assisted by Google Antigravity AI</i>", affiliation_style))
    
    # Abstract
    story.append(Paragraph("ABSTRACT", abstract_heading))
    abstract_text = (
        "High-velocity kinematic tails in cold stellar streams are widely sought as gravitational signatures of dark matter subhalo perturbations. "
        "Recent analyses of DESI DR2 spectroscopy in the GD-1 stream reported an apparent kinematic tail excess. Here, we present an end-to-end "
        "comprehensive kinematic audit across official DESI DR2 data products (pre-release Table 7 and the updated multi-component catalog, "
        "Jarvis et al. 2026, arXiv:2604.20958; Zenodo concept DOI: 10.5281/zenodo.19258455). First, conventional hard probability thresholding (P_thin &gt;= 0.90) "
        "reveals severe sensitivity to catalog classification choices: a single borderline star (Star 28) shifting by delta P = -0.0166 flips the inferred significance "
        "by 12.5x across alpha = 0.05 (p = 0.0043 to 0.0539). Second, replacing step functions with a threshold-free continuous Bayesian score (T_cont) demonstrates "
        "that continuous tail statistics remain detectable (T_obs = 27.82, p = 0.0012) and persist even after omitting Star 28 (p = 0.0113); an audit of Star 28's 14 immediate neighbors (+/- 3 deg) reveals they are dynamically cold "
        "(-1.40 km/s, -0.25 sigma), showing no evidence for a localized single-subhalo collision dipole. Third, cross-matching DESI DR2 kinematics with Gaia DR3 proper motions "
        "reveals extreme velocity anisotropy: the line-of-sight velocity variance is 6.2x more excited than 2D tangential proper motions (R_aniso = 0.1613, r = 0.0464), "
        "directly contradicting 3D isotropic heating. Exactly 9 stars (78.5% of tail energy) possess cold proper motions (z_tan &lt;= 1.36 sigma), nominal "
        "single-star Gaia astrometry (RUWE ~ 1.01), and stream metallicity ([Fe/H] &lt;= -1.81), causing the residual tail score to collapse to null expectation (p_residual = 0.7727). "
        "Fourth, archival SDSS DR17 spectra show multi-epoch velocity variability across surveys: thin stream member Gaia DR3 804386042407893376 in the spur bin exhibited "
        "delta V = +28.64 +/- 7.17 km/s (~3.5 sigma calibrated against empirical cross-survey dispersion, an independent recovery of the ~30 km/s variable noted in Section II of Nibauer et al. 2026), "
        "while Outlier 13 in the diffuse cocoon envelope showed delta V = +20.60 +/- 5.88 km/s (~2.4 sigma calibrated), Star 4 showed a modest drift (-5.79 +/- 2.44 km/s, ~1.4 sigma calibrated, consistent with noise), "
        "and Star 28 is dynamically constant (0.24 sigma). Fifth, direct bin-by-bin confrontation with the multi-survey modeling of Nibauer et al. (2026) reveals that DESI DR2 "
        "is consistent with their finding that the gap/spur bin ([-40, -20] deg) is the hottest segment along the stream (3.32 +/- 0.49 km/s raw, reducing to 2.72 +/- 0.41 km/s, or 2.68 +/- 0.42 km/s without the spur candidate), "
        "though the elevation over the trailing (1.87 +/- 0.32 km/s) and leading (1.57 +/- 0.30 km/s) arms is modest (~1.5 - 2 sigma). "
        "We conclude that single-epoch velocity tails in GD-1 are heavily contaminated by binaries and envelope stars, while the localized spur excess could reflect subhalo perturbations, "
        "cocoon contamination, or undetected binaries (exemplified by the 29 km/s velocity shift of Gaia 804386042407893376). Dedicated dynamical forward-modeling is required before robust "
        "fundamental physics bounds can be inferred."
    )
    story.append(Paragraph(abstract_text, abstract_style))
    story.append(Spacer(1, 2))
    
    # Section 1
    story.append(Paragraph("1. Introduction", h1_style))
    story.append(Paragraph(
        "Cold stellar streams in the Galactic halo serve as sensitive gravitational balances for dark matter subhalos. Encounters with intervening substructure "
        "are expected to impart localized velocity kicks and density gaps (Carlberg 2012; Erkal &amp; Belokurov 2015; Bonaca et al. 2019). Recently, Jarvis et al. (2026; "
        "arXiv:2604.20958; Zenodo concept DOI: 10.5281/zenodo.19258455) modeled GD-1 DESI DR2 kinematics, while Carlberg et al. (2026, arXiv:2606.15863) modeled joint "
        "stream dynamics. Stellar binaries in streams present a classic contaminant that inflates apparent dispersions; recently, Nibauer et al. (2026, ApJ, 1004, 62) "
        "used repeat spectroscopy across four surveys (DESI EDR, SDSS DR9, LAMOST DR8, and MMT Hectochelle) to model the GD-1 binary population and constrain subhalo heating. "
        "Notably, Nibauer et al. already incorporated early DESI data (EDR); our work provides a kinematic audit of the complete DESI DR2 GD-1 sample, examining probability "
        "threshold fragility, continuous Bayesian scoring, transverse proper-motion orthogonality, empirical SDSS-DESI zero-points, and direct bin-by-bin confrontation with Nibauer et al.",
        body_style
    ))
    
    # Section 2
    story.append(Paragraph("2. Probability Threshold Sensitivity and Star 28 (v43)", h1_style))
    story.append(Paragraph(
        "In the pre-release catalog (v3), filtering stars by P_thin &gt;= 0.90 yielded 306 stream stars, including Star 28 (Gaia DR3 804436551223805568, phi1 = -32.93 deg, delta V_GSR = -17.88 km/s, "
        "P_thin = 0.9106), giving a tail excess p = 0.0043. In the updated catalog (v4), truncation re-normalization adjusted Star 28 to P_thin = 0.8940 (-1.8% change). "
        "Applying the step function deleted Star 28, multiplying the p-value by 12.5x to p = 0.0539. Omitting Star 28 collapsed the p-value to p = 0.3980.",
        body_style
    ))
    story.append(Paragraph(
        "Furthermore, at P_thin &gt;= 0.80, the observed count was N_obs = 316, strictly below the baseline 95% predictive interval [321, 330] (p_N = 0.0006), "
        "revealing that fixed-fraction binomial forward models failed spatial density gradients along the stream.",
        body_style
    ))
    
    # Section 3
    story.append(Paragraph("3. Continuous Bayesian Scoring &amp; Spatial Audit (v44)", h1_style))
    story.append(Paragraph(
        "To eliminate step-function artifacts, we introduced the continuous score T_cont = sum(P_thin_i * max(0, z_i^2 - 2.5^2)). "
        "Under T_cont, Star 28 adjusts smoothly from 5.84 to 5.72. Across 10,000 null realizations with rejection truncation, T_obs = 27.8212 yields p = 0.0012, "
        "and excluding Star 28 yields p = 0.0113. The velocity tail is statistically robust. Star 28 is dynamically isolated: its 14 immediate neighbors within +/- 3 deg "
        "are dynamically cold (-1.40 km/s, -0.25 sigma), showing no evidence for a localized single-subhalo collision dipole. Field-level Dirichlet-Multinomial overdispersion "
        "(sigma_pi = 0.1536) widened the 95% predictive interval at P_thin &gt;= 0.80 to [309, 329], reconciling N_obs = 316 (p_N = 0.5400).",
        body_style
    ))
    
    story.append(PageBreak())
    
    # -------------------------------------------------------------------------
    # PAGE 2: Section 4, Figure 1 (v45), Table 1
    # -------------------------------------------------------------------------
    story.append(Paragraph("4. Physical Attribution &amp; Velocity Anisotropy (v45)", h1_style))
    story.append(Paragraph(
        "Isotropic 3D dark-matter subhalo heating predicts R_aniso = &lt;z_tan^2 / 2&gt; / &lt;z_v^2&gt; ~ 1.0. Empirical evaluation of the 14 tail stars reveals "
        "&lt;z_v^2&gt; = 11.2112, &lt;z_tan^2 / 2&gt; = 1.8085, giving R_aniso = 0.1613. Radial velocity is 6.2x more excited than proper motions, with r = 0.0464 (p = 0.8749). "
        "A KS test on proper motions against the cold core yields p_KS = 0.2485.",
        body_style
    ))
    story.append(Paragraph(
        "Exactly 9 stars (78.5% of tail energy) are candidate binaries: cold proper motions (z_tan &lt;= 1.36 sigma), nominal RUWE (~1.012), and GD-1 metallicity "
        "([Fe/H] &lt;= -1.81). At the 8-10 kpc distance of GD-1, an orbital separation of 1-2 AU subtends only ~0.1-0.2 mas, well below Gaia photocentric wobble sensitivity; moreover, equal-luminosity pairs exhibit no photocenter shift, and wide binaries show negligible acceleration over Gaia's baseline. Nominal RUWE thus neither supports nor rules out binarity at this distance; the decisive physical discriminant is radial velocity anisotropy. Subtracting binaries leaves residual score T_res = 5.9901. Against the null distribution (9.67 +/- 4.70), p_residual = 0.7727.",
        body_style
    ))
    
    # Figure 1 (v45)
    story.append(Spacer(1, 4))
    if os.path.exists("D:/ARGON/v45_physical_attribution.png"):
        story.append(Image("D:/ARGON/v45_physical_attribution.png", width=504, height=173.6))
        story.append(Paragraph("<b>Figure 1:</b> Physical Attribution Audit (v45). (a) 2D kinematic orthogonality (z_tan vs |z_v|). (b) Gaia DR3 astrometric RUWE (~1.01). (c) DESI metallicity along phi1. (d) Tail score energy decomposition showing collapse to p = 0.7727.", caption_style))
    
    # Table of 14 stars
    story.append(Spacer(1, 4))
    story.append(Paragraph("Table 1: The 14 Active Tail Stars in Zenodo v4 and Multi-Epoch Spectroscopic Shifts", h2_style))
    table_data = [
        [Paragraph(h, table_header) for h in ["#", "Source ID", "phi1", "dV_GSR", "z_v", "z_tan", "RUWE", "[Fe/H]", "T_i", "Type", "Multi-Epoch Shift"]]
    ]
    
    raw_table = [
        ["1", "804436551223805568 (Star 28)", "-32.93", "-17.88", "-3.56", "0.65", "0.996", "-2.37", "5.72", "A", "+2.40 +/- 8.43 km/s (0.29s)"],
        ["2", "637140290779827072", "-57.93", "-18.85", "-4.30", "1.36", "0.910", "-2.64", "3.79", "A", "---"],
        ["3", "1577287126973726592", "-3.25", "-25.39", "-3.49", "1.11", "1.038", "-2.21", "3.63", "A", "---"],
        ["4", "794510439424008832", "-42.26", "+8.31", "+3.12", "1.08", "1.069", "-2.67", "2.76", "A", "---"],
        ["5", "831423209431680384", "-24.55", "-18.30", "-3.63", "2.76", "1.007", "-2.39", "2.26", "B", "---"],
        ["6", "829645264771612800 (Star 4)", "-27.35", "-7.85", "-2.89", "1.16", "1.098", "-2.03", "1.94", "A", "-5.79 +/- 2.44 km/s (2.37s, ~1.4s emp)"],
        ["7", "1566540564752850816", "-1.90", "+22.88", "+3.78", "2.51", "0.984", "-2.29", "1.90", "B", "---"],
        ["8", "794347947923184640", "-43.79", "-14.52", "-4.05", "1.05", "1.032", "-2.12", "1.64", "A", "---"],
        ["9", "803361336230087936", "-36.14", "-10.39", "-3.57", "2.05", "0.988", "-2.66", "1.58", "B", "+2.06 +/- 9.16 km/s (0.22s)"],
        ["10", "839367829479258240", "-18.48", "+10.55", "+2.84", "0.44", "1.000", "-1.81", "1.16", "A", "---"],
        ["11", "1577024756011421696", "-2.74", "+10.03", "+2.74", "0.88", "0.965", "-1.96", "1.07", "A", "---"],
        ["12", "803407064742685568", "-36.24", "+9.67", "+2.77", "2.50", "1.048", "-2.21", "0.17", "B", "---"],
        ["13", "796893149840353408 (Outlier 13)", "-39.46", "+6.80", "+2.57", "1.02", "1.009", "-2.68", "0.12", "A", "+20.60 +/- 5.88 km/s (3.50s, ~2.4s emp)"],
        ["14", "840450504835587328", "-14.60", "-20.08", "-2.99", "4.14", "1.020", "-2.54", "0.08", "B", "---"]
    ]
    for row in raw_table:
        table_data.append([Paragraph(cell, table_cell) for cell in row])
        
    t_obj = Table(table_data, colWidths=[16, 114, 34, 44, 32, 32, 38, 34, 26, 30, 104])
    t_obj.setStyle(TableStyle([
        ('BACKGROUND', (0,0), (-1,0), colors.HexColor("#1A365D")),
        ('ALIGN', (0,0), (-1,-1), 'CENTER'),
        ('VALIGN', (0,0), (-1,-1), 'MIDDLE'),
        ('GRID', (0,0), (-1,-1), 0.4, colors.HexColor("#CBD5E0")),
        ('ROWBACKGROUNDS', (0,1), (-1,-1), [colors.white, colors.HexColor("#F7FAFC")]),
        ('TOPPADDING', (0,0), (-1,-1), 1.2),
        ('BOTTOMPADDING', (0,0), (-1,-1), 1.2),
    ]))
    story.append(t_obj)
    
    story.append(PageBreak())
    
    # -------------------------------------------------------------------------
    # PAGE 3: Section 5, Figure 2 (v46)
    # -------------------------------------------------------------------------
    story.append(Paragraph("5. Multi-Epoch Spectroscopic Audit &amp; Radial-Velocity Variability (v46)", h1_style))
    story.append(Paragraph(
        "Cross-matching against SDSS DR17 spectroscopy identified independent multi-epoch observations for evaluating radial-velocity variability:",
        body_style
    ))
    story.append(Paragraph(
        "<b>Outlier 13 (Gaia DR3 796893149840353408):</b> Classified predominantly as a cocoon member (P_cocoon = 0.67, P_thin = 0.33) at phi1 = -39.46 deg. "
        "In SDSS DR17 (2016, MJD 57450, Plate 8854, Fiber 508, specObjID 9968857538152912896), V_LOS = -32.60 +/- 5.76 km/s. In DESI DR2 (2021), "
        "V_LOS = -12.00 +/- 1.18 km/s. The raw shift is delta V = +20.60 +/- 5.88 km/s (3.50 sigma formal). Calibrating against empirical cross-survey dispersion "
        "(~7.0 km/s) and zero-point (+4.1 km/s) yields a net shift of +16.54 km/s (~2.37 sigma), consistent with a candidate binary in the cocoon envelope.",
        bullet_style
    ))
    story.append(Paragraph(
        "<b>Star 4 (Gaia DR3 829645264771612800):</b> At phi1 = -27.35 deg, Star 4 sat inside the cold core in SDSS DR17 (2009, MJD 54884, Plate 3258, Fiber 371, "
        "specObjID 3668283958136774656) at V_LOS = -100.26 +/- 2.31 km/s. In DESI DR2 (2021), it was measured at -106.05 +/- 0.79 km/s. This yields a raw drift of "
        "delta V = -5.79 +/- 2.44 km/s (2.37 sigma formal, ~1.41 sigma calibrated).",
        bullet_style
    ))
    story.append(Paragraph(
        "<b>Cocoon Star 9 (Gaia DR3 803361336230087936):</b> In SDSS DR17 (2004, MJD 53033, Plate 1356, Fiber 106, specObjID 1526749461622056960), V_LOS = -52.78 +/- 9.06 km/s "
        "vs DESI -50.72 +/- 1.33 km/s (delta V = +2.06 +/- 9.16 km/s, 0.22 sigma). Its velocity shows no detectable variability over 18 years, consistent with a non-variable member of the diffuse cocoon.",
        bullet_style
    ))
    story.append(Paragraph(
        "<b>Star 28 (Gaia DR3 804436551223805568):</b> In SDSS DR17 (2009, MJD 54941, Plate 3287, Fiber 94, specObjID 3700858915211288576), Star 28 was measured at "
        "V_LOS = -81.10 +/- 7.27 km/s vs DESI -78.70 +/- 4.27 km/s (delta V = +2.40 +/- 8.43 km/s, 0.29 sigma formal, 0.24 sigma calibrated), exhibiting no detectable variation over 12 years.",
        bullet_style
    ))
    story.append(Paragraph(
        "<b>Systematic Zero-Points &amp; Empirical Calibration:</b> Direct cross-matching of DESI DR2 stars against SDSS DR17 on SkyServer yields 122 matches; for 68 clean "
        "non-binary stars (|delta V| &lt; 15 km/s), the median offset is +4.06 km/s (mean +3.51 +/- 0.85 km/s, std 6.98 km/s). This ~7.0 km/s empirical scatter implies formal "
        "cross-survey errors are underestimated. When evaluated against this scatter, Outlier 13 is ~2.37 sigma and Star 4 is ~1.41 sigma. Crucially, thin stream member "
        "Gaia DR3 804386042407893376 in the spur bin (phi1 = -33.77 deg, P_thin = 0.961) shows delta V = +28.64 +/- 7.17 km/s between SDSS (MJD 55570, Plate 4562, Fiber 865, "
        "specObjID 5136593237855983616) and DESI. Calibrated against empirical scatter, this shift is ~3.52 sigma, making it the strongest variable candidate and possibly the same "
        "star noted by Nibauer et al. (2026). This star was not in our 14 tail stars and was not removed from the spur bin, demonstrating that active binaries remain embedded within the clean core.",
        bullet_style
    ))
    
    # Figure 2 (v46)
    story.append(Spacer(1, 6))
    if os.path.exists("D:/ARGON/v46_orbital_phase_confirmation.png"):
        story.append(Image("D:/ARGON/v46_orbital_phase_confirmation.png", width=504, height=195.0))
        story.append(Paragraph("<b>Figure 2:</b> Multi-Epoch Spectroscopic Check (v46). (a) Outlier 13 zero-point crossover (+20.6 km/s, 3.50 sigma) and Star 4 cold-core exit (-5.8 km/s, 2.37 sigma). (b) Epistemic resolution progression across v43 to v46.", caption_style))

    story.append(PageBreak())

    # -------------------------------------------------------------------------
    # PAGE 4: Section 6, Figure 3 (v47)
    # -------------------------------------------------------------------------
    story.append(Paragraph("6. Spatial Kinematics and Confrontation with Stream Perturbations (v47)", h1_style))
    story.append(Paragraph(
        "Subtracting the 9 candidate Type A binaries reveals a dynamically cold stream core (N = 536), "
        "enabling a robust evaluation of intrinsic dispersion and direct spatial confrontation with published stream heating models:",
        body_style
    ))
    story.append(Paragraph(
        "<b>Global Stream Core Velocity Dispersion:</b> Evaluated across the entire [-60, 0] deg window (including the spur), the weighted Bayesian profile likelihood "
        "yields sigma_int = 1.88 (+0.19 / -0.18) km/s, with a 95% credible upper bound of sigma_int,95% &lt;= 2.27 km/s (a 17.1% like-for-like reduction relative to the raw profile likelihood of 2.27 km/s, or 18.9% under continuous weights; this reduction is an upper bound because candidates were selected partly on velocity). "
        "The spatial breakdown across the three bins shows that this global value is an intermediate compromise between the colder outer arms (1.57 - 1.87 km/s) and the spur bin (2.72 km/s).",
        bullet_style
    ))
    story.append(Paragraph(
        "<b>Bin-by-Bin Confrontation with Nibauer et al. (2026):</b> Evaluating DESI DR2 across the three Nibauer et al. (2026) bins reveals: "
        "Trailing ([-60, -40] deg): 2.33 +/- 0.39 to 1.87 +/- 0.32 km/s (-19.6%); Gap/Spur ([-40, -20] deg): 3.32 +/- 0.49 to 2.72 +/- 0.41 km/s (-17.9%, or 2.68 +/- 0.42 km/s without spur variable Gaia 804386042407893376); "
        "Leading ([-20, 0] deg): 1.73 +/- 0.31 to 1.57 +/- 0.30 km/s (-9.1%). While the spur remains the hottest bin along the stream, the difference over the trailing (1.87 km/s) and leading (1.57 km/s) arms is modest (~1.5 - 2 sigma). "
        "Possible contributors to the gap between our cleaned spur dispersion (2.72 km/s) and Nibauer et al.'s ~4.8 km/s include: (i) Sample selection &amp; weighting: Nibauer et al. adopt an unweighted |delta v| &lt; 30 km/s window. Applying this exact unweighted cut to our DESI DR2 spur sample yields sigma = 6.49 +/- 0.71 km/s, reflecting substantial envelope contamination; conversely, because P_thin incorporates velocity, weighting by P_thin pulls the dispersion down by construction (3.32 km/s raw, 2.72 km/s cleaned). "
        "(ii) Binary treatment: their single parametric jitter term sigma_binary treats binaries globally, whereas our audit identifies individual candidates (see Phillips et al. 2026 for stream binary populations). We emphasize we have not reproduced their measurement on a common sample.",
        bullet_style
    ))
    story.append(Paragraph(
        "<b>Origin of the Spur Excess and Future Modeling:</b> The residual kinematic excess in the spur bin (2.72 +/- 0.41 km/s, or 2.68 w/o 804) admits three physical interpretations: "
        "(i) localized subhalo heating (Nibauer et al. 2026), (ii) diffuse cocoon envelope stars (Jarvis et al. 2026), or (iii) undetected spectroscopic binaries. "
        "The presence of hidden binaries is directly confirmed by Gaia DR3 804386042407893376, an independent recovery of the ~30 km/s variable noted in Section II of Nibauer et al. (2026). "
        "In DESI DR2, its track residual (+7.47 km/s) sits within 1.5 sigma of the combined track envelope (sqrt(4.24^2 + 2.7^2) = 5.0 km/s), explaining why it was retained in the clean core. "
        "Disentangling these mechanisms will require dedicated multi-epoch radial-velocity monitoring and tailored forward dynamical simulations.",
        bullet_style
    ))
    
    # Figure 3 (v47)
    story.append(Spacer(1, 4))
    if os.path.exists("D:/ARGON/v47_kinematic_confrontation.png"):
        story.append(Image("D:/ARGON/v47_kinematic_confrontation.png", width=504, height=177.5))
        story.append(Paragraph("<b>Figure 3:</b> Kinematic Profile and Spatial Confrontation of the GD-1 Stream in DESI DR2 (v47). (a) Weighted Bayesian profile likelihood for sigma_int (1.88 km/s, 17.1% reduction vs raw profile fit of 2.27 km/s; upper bound because candidates were selected partly on velocity). (b) Spatial kinematic confrontation across the three Nibauer et al. (2026) bins, showing the spur bin drop (3.32 to 2.72 km/s, and 2.68 km/s w/o 804) while outer arms remain cold (1.57 - 1.87 km/s). (c) Velocity residuals along phi1. (d) Multi-epoch velocity shifts across SDSS DR17 and DESI DR2, highlighting independent recovery of spur candidate binary Gaia DR3 804386042407893376 (delta V = +28.64 km/s, ~3.5 sigma_emp), Outlier 13 (~2.4 sigma_emp), and Star 4 (1.4 sigma_emp, consistent with noise).", caption_style))

    story.append(PageBreak())

    # -------------------------------------------------------------------------
    # PAGE 5: Section 7 (Conclusions), Acknowledgments, References
    # -------------------------------------------------------------------------
    story.append(Paragraph("7. Discussion and Conclusions", h1_style))
    story.append(Paragraph(
        "The end-to-end scientific progression across v43 through v47 establishes four core methodological and astrophysical principles:",
        body_style
    ))
    story.append(Paragraph(
        "1. <b>Continuous Scoring Protects Against Boundary Cliff Artifacts:</b> Step-function thresholding (P &gt;= 0.90) creates severe fragility. Continuous scoring (T_cont) provides monotonic, stable inference.",
        bullet_style
    ))
    story.append(Paragraph(
        "2. <b>Proper Motions Discriminate Binaries from Dark Matter Heating:</b> Subhalos heat streams isotropically (R_aniso ~ 1.0). Extreme line-of-sight anisotropy (R_aniso = 0.1613) combined with cold proper motions (z_tan &lt;= 1.36 sigma) is consistent with binaries.",
        bullet_style
    ))
    story.append(Paragraph(
        "3. <b>Multi-Epoch Spectroscopy Identifies Radial-Velocity Shifts:</b> Single-epoch spectroscopic velocity tails cannot be attributed to dark matter perturbations without multi-epoch verification. Calibrated against empirical cross-survey dispersion (~7.0 km/s), thin stream member Gaia DR3 804386042407893376 in the spur bin exhibited a delta V = +28.64 +/- 7.17 km/s velocity shift (~3.52 sigma, an independent recovery of the ~30 km/s variable noted in Section II of Nibauer et al. 2026), while Outlier 13 in the cocoon showed ~2.37 sigma, Star 4 showed modest drift (~1.41 sigma), and Star 28 is dynamically constant (0.24 sigma).",
        bullet_style
    ))
    story.append(Paragraph(
        "4. <b>Confrontation with Stream Heating Models and Origin of Spur Excess:</b> Subtracting candidate binaries reduces GD-1's global core velocity dispersion to sigma_int = 1.88 (+0.19 / -0.18) km/s (95% upper limit 2.26 km/s). Direct confrontation with the multi-survey bins of Nibauer et al. (2026) reveals that DESI DR2 is consistent with their peak dispersion in the gap/spur bin ([-40, -20] deg, 3.32 to 2.72 km/s), the spur remains the hottest bin along the stream, though the elevation over the outer arms (1.57 - 1.87 km/s) is modest (~1.5 - 2 sigma). This localized spur excess could stem from compact subhalo perturbations, diffuse cocoon envelope contamination, or undetected binaries (exemplified by the 29 km/s velocity shift of Gaia 804386042407893376). Dedicated dynamical forward-modeling is required before robust fundamental physics bounds can be inferred.",
        bullet_style
    ))
    
    story.append(Spacer(1, 2))
    story.append(Paragraph("Acknowledgments", h2_style))
    story.append(Paragraph(
        "This research made use of public data products from the DESI DR2 survey, the European Space Agency Gaia mission, "
        "and the Sloan Digital Sky Survey (SDSS) DR17. Autonomous kinematic workflow orchestration, Bayesian profile "
        "likelihood calculations, and data visualization were conducted using Google DeepMind's Antigravity AI research agent.",
        body_style
    ))
    
    story.append(Spacer(1, 2))
    story.append(Paragraph("References", h2_style))
    refs = [
        "Bonaca, A., Hogg, D. W., Price-Whelan, A. M., et al. 2019, ApJ, 880, 38",
        "Bovy, J., Erkal, D., &amp; Sanders, J. L. 2017, MNRAS, 466, 628",
        "Carlberg, R. G. 2012, ApJ, 748, 20",
        "Carlberg, R. G., Li, T. S., Jarvis, E., et al. 2026, arXiv:2606.15863",
        "Cooper, A. P., Koposov, S. E., Dey, A., et al. 2023, ApJ, 947, 37",
        "Erkal, D. &amp; Belokurov, V. 2015, MNRAS, 450, 1136",
        "Garrison-Kimmel, S., Wetzel, A., Bullock, J. S., et al. 2017, MNRAS, 471, 1709",
        "Grillmair, C. J. &amp; Dionatos, O. 2006, ApJL, 643, L17",
        "Helmi, A. &amp; White, S. D. M. 2001, MNRAS, 307, 495",
        "Jarvis, E., Li, T. S., Koposov, S. E., et al. 2026, arXiv:2604.20958, Zenodo, doi:10.5281/zenodo.19258455",
        "Johnston, K. V., Hernquist, L., &amp; Bolte, M. 1996, ApJ, 465, 278",
        "Kelley, T., Bullock, J. S., Garrison-Kimmel, S., et al. 2019, MNRAS, 487, 4409",
        "Koposov, S. E., Rix, H.-W., &amp; Hogg, D. W. 2010, ApJ, 712, 260",
        "Koposov, S. E., Cooper, A. P., Dey, A., et al. 2024, MNRAS, 528, 6934",
        "Malhan, K., Ibata, R. A., &amp; Martin, N. F. 2018, MNRAS, 481, 3442",
        "Nibauer, J., Bonaca, A., Price-Whelan, A. M., Spergel, D. N., &amp; Greene, J. E. 2026, ApJ, 1004, 62",
        "Phillips, N., Conroy, C., Nibauer, J., et al. 2026, arXiv:2603.06790",
        "Price-Whelan, A. M. &amp; Bonaca, A. 2018, ApJL, 863, L20",
        "Schneider, A., Smith, R. E., Maccio, A. V., &amp; Moore, B. 2012, MNRAS, 424, 684",
        "Thomas, G. F., Famaey, B., Ibata, R., et al. 2018, A&amp;A, 609, A44"
    ]
    for r in refs:
        story.append(Paragraph(r, ParagraphStyle('RefStyle', parent=styles['Normal'], fontName='Helvetica', fontSize=6.6, leading=8.1, textColor=colors.HexColor("#4A5568"), leftIndent=10, spaceAfter=1.2)))

    doc.build(story, canvasmaker=NumberedCanvas)
    print(f"Successfully compiled PDF: {PDF_FILENAME}")

if __name__ == '__main__':
    build_pdf()
