#!/usr/bin/env python3
"""Build the oscillon research orientation PDF from report-source.md."""

from __future__ import annotations

import html
import re
from pathlib import Path

from reportlab.graphics.shapes import Circle, Drawing, Line, Polygon, Rect, String
from reportlab.lib import colors
from reportlab.lib.colors import HexColor
from reportlab.lib.enums import TA_CENTER, TA_LEFT
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
from reportlab.lib.units import mm
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.platypus import (
    BaseDocTemplate,
    Flowable,
    Frame,
    HRFlowable,
    KeepTogether,
    ListFlowable,
    ListItem,
    PageBreak,
    PageTemplate,
    Paragraph,
    Spacer,
    Table,
    TableStyle,
)
from reportlab.platypus.tableofcontents import TableOfContents


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "report-source.md"
OUTPUT = ROOT / "output" / "pdf" / "oscillons_research_orientation.pdf"

NAVY = HexColor("#10243E")
BLUE = HexColor("#2F6690")
TEAL = HexColor("#2A9D8F")
GOLD = HexColor("#E9C46A")
CORAL = HexColor("#E76F51")
INK = HexColor("#17212B")
MUTED = HexColor("#516172")
PALE_BLUE = HexColor("#EEF5F9")
PALE_TEAL = HexColor("#EAF6F3")
PALE_GOLD = HexColor("#FBF5E4")
PALE_CORAL = HexColor("#FCEFEA")
RULE = HexColor("#CBD6DF")


def register_fonts() -> None:
    font_files = {
        "Body": "/usr/share/fonts/google-crosextra-caladea-fonts/Caladea-Regular.ttf",
        "Body-Bold": "/usr/share/fonts/google-crosextra-caladea-fonts/Caladea-Bold.ttf",
        "Body-Italic": "/usr/share/fonts/google-crosextra-caladea-fonts/Caladea-Italic.ttf",
        "Body-BoldItalic": "/usr/share/fonts/google-crosextra-caladea-fonts/Caladea-BoldItalic.ttf",
        "Sans": "/usr/share/fonts/google-carlito-fonts/Carlito-Regular.ttf",
        "Sans-Bold": "/usr/share/fonts/google-carlito-fonts/Carlito-Bold.ttf",
        "Sans-Italic": "/usr/share/fonts/google-carlito-fonts/Carlito-Italic.ttf",
        "Sans-BoldItalic": "/usr/share/fonts/google-carlito-fonts/Carlito-BoldItalic.ttf",
        "Mono": "/usr/share/fonts/adwaita-mono-fonts/AdwaitaMono-Regular.ttf",
        "Mono-Bold": "/usr/share/fonts/adwaita-mono-fonts/AdwaitaMono-Bold.ttf",
    }
    for name, path in font_files.items():
        pdfmetrics.registerFont(TTFont(name, path))
    pdfmetrics.registerFontFamily(
        "Body", normal="Body", bold="Body-Bold", italic="Body-Italic", boldItalic="Body-BoldItalic"
    )
    pdfmetrics.registerFontFamily(
        "Sans", normal="Sans", bold="Sans-Bold", italic="Sans-Italic", boldItalic="Sans-BoldItalic"
    )


def make_styles():
    sample = getSampleStyleSheet()
    styles = {}
    styles["Body"] = ParagraphStyle(
        "Body",
        parent=sample["BodyText"],
        fontName="Body",
        fontSize=9.35,
        leading=12.25,
        textColor=INK,
        spaceAfter=5.8,
        splitLongWords=True,
    )
    styles["Small"] = ParagraphStyle(
        "Small",
        parent=styles["Body"],
        fontName="Sans",
        fontSize=7.7,
        leading=9.5,
        textColor=MUTED,
        spaceAfter=4,
    )
    styles["Heading1"] = ParagraphStyle(
        "Heading1",
        parent=sample["Heading1"],
        fontName="Sans-Bold",
        fontSize=18,
        leading=21.5,
        textColor=NAVY,
        spaceBefore=13,
        spaceAfter=7,
        keepWithNext=True,
    )
    styles["Heading2"] = ParagraphStyle(
        "Heading2",
        parent=sample["Heading2"],
        fontName="Sans-Bold",
        fontSize=12.3,
        leading=14.5,
        textColor=BLUE,
        spaceBefore=9,
        spaceAfter=4.5,
        keepWithNext=True,
    )
    styles["Heading3"] = ParagraphStyle(
        "Heading3",
        parent=sample["Heading3"],
        fontName="Sans-Bold",
        fontSize=10.1,
        leading=12.2,
        textColor=TEAL,
        spaceBefore=7,
        spaceAfter=3,
        keepWithNext=True,
    )
    styles["Executive"] = ParagraphStyle(
        "Executive",
        parent=styles["Heading1"],
        fontSize=20,
        leading=23,
        textColor=NAVY,
        spaceBefore=0,
        spaceAfter=9,
    )
    styles["Callout"] = ParagraphStyle(
        "Callout",
        parent=styles["Body"],
        fontName="Sans",
        fontSize=9.25,
        leading=12.3,
        spaceAfter=0,
    )
    styles["Equation"] = ParagraphStyle(
        "Equation",
        parent=styles["Body"],
        fontName="Mono",
        fontSize=8.25,
        leading=11.2,
        alignment=TA_CENTER,
        textColor=NAVY,
        spaceAfter=0,
    )
    styles["Table"] = ParagraphStyle(
        "Table",
        parent=styles["Body"],
        fontName="Sans",
        fontSize=7.45,
        leading=9.1,
        spaceAfter=0,
    )
    styles["TableHead"] = ParagraphStyle(
        "TableHead",
        parent=styles["Table"],
        fontName="Sans-Bold",
        textColor=colors.white,
        leading=9.3,
    )
    styles["List"] = ParagraphStyle(
        "List",
        parent=styles["Body"],
        leftIndent=0,
        firstLineIndent=0,
        spaceAfter=2.2,
    )
    styles["Caption"] = ParagraphStyle(
        "Caption",
        parent=styles["Small"],
        fontName="Sans-Italic",
        alignment=TA_CENTER,
        textColor=MUTED,
        spaceBefore=3,
        spaceAfter=7,
    )
    styles["CoverKicker"] = ParagraphStyle(
        "CoverKicker",
        parent=sample["BodyText"],
        fontName="Sans-Bold",
        fontSize=9,
        leading=11,
        textColor=GOLD,
        spaceAfter=10,
        uppercase=True,
    )
    styles["CoverTitle"] = ParagraphStyle(
        "CoverTitle",
        parent=sample["Title"],
        fontName="Sans-Bold",
        fontSize=29,
        leading=31.5,
        textColor=colors.white,
        alignment=TA_LEFT,
        spaceAfter=12,
    )
    styles["CoverSub"] = ParagraphStyle(
        "CoverSub",
        parent=sample["BodyText"],
        fontName="Body",
        fontSize=13.5,
        leading=18,
        textColor=HexColor("#DBE7F0"),
        spaceAfter=16,
    )
    styles["CoverMeta"] = ParagraphStyle(
        "CoverMeta",
        parent=sample["BodyText"],
        fontName="Sans",
        fontSize=9.2,
        leading=13,
        textColor=HexColor("#B8CBDA"),
    )
    return styles


def inline_markup(text: str) -> str:
    """Convert a restrained subset of Markdown to ReportLab paragraph markup."""
    tokens: list[str] = []

    def token(value: str) -> str:
        idx = len(tokens)
        tokens.append(value)
        return f"@@TOKEN{idx}@@"

    def link_repl(match: re.Match) -> str:
        label = html.escape(match.group(1), quote=False)
        url = html.escape(match.group(2), quote=True)
        return token(f'<a href="{url}" color="#2F6690"><u>{label}</u></a>')

    def code_repl(match: re.Match) -> str:
        code = html.escape(match.group(1), quote=False)
        return token(f'<font name="Mono" color="#244A64">{code}</font>')

    text = re.sub(r"\[([^\]]+)\]\(([^)]+)\)", link_repl, text)
    text = re.sub(r"`([^`]+)`", code_repl, text)
    text = html.escape(text, quote=False)
    text = re.sub(r"\*\*(.+?)\*\*", r"<b>\1</b>", text)
    text = re.sub(r"(?<!\*)\*([^*]+?)\*(?!\*)", r"<i>\1</i>", text)
    for idx, value in enumerate(tokens):
        text = text.replace(f"@@TOKEN{idx}@@", value)
    return text


def paragraph(text: str, style, bullet: str | None = None) -> Paragraph:
    return Paragraph(inline_markup(text), style, bulletText=bullet)


class ResearchDocTemplate(BaseDocTemplate):
    def __init__(self, filename: str, styles, **kwargs):
        super().__init__(filename, **kwargs)
        self.styles = styles
        self._heading_counter = 0

    def afterFlowable(self, flowable: Flowable) -> None:
        if not isinstance(flowable, Paragraph):
            return
        style_name = flowable.style.name
        if style_name not in {"Heading1", "Heading2"}:
            return
        level = 0 if style_name == "Heading1" else 1
        text = flowable.getPlainText()
        key = f"heading-{self._heading_counter}"
        self._heading_counter += 1
        self.canv.bookmarkPage(key)
        self.canv.addOutlineEntry(text, key, level=level, closed=False)
        self.notify("TOCEntry", (level, text, self.page, key))


def cover_background(canvas, doc) -> None:
    width, height = A4
    page = canvas.getPageNumber()
    canvas.saveState()
    if page == 1:
        canvas.setFillColor(NAVY)
        canvas.rect(0, 0, width, height, stroke=0, fill=1)
        canvas.setFillColor(HexColor("#173A59"))
        canvas.circle(width * 0.82, height * 0.83, 78 * mm, stroke=0, fill=1)
        canvas.setFillColor(HexColor("#1C4967"))
        canvas.circle(width * 0.82, height * 0.83, 55 * mm, stroke=0, fill=1)
        canvas.setStrokeColor(TEAL)
        canvas.setLineWidth(2.2)
        points = []
        for n in range(130):
            x = 16 * mm + n * (178 * mm / 129)
            envelope = 23 * mm * (2.71828 ** (-((n - 74) / 31) ** 2))
            y = 62 * mm + envelope * __import__("math").sin(n * 0.43)
            points.extend([x, y])
        canvas.line(16 * mm, 62 * mm, 194 * mm, 62 * mm)
        for i in range(0, len(points) - 2, 2):
            canvas.line(points[i], points[i + 1], points[i + 2], points[i + 3])
        canvas.setFillColor(GOLD)
        canvas.circle(161 * mm, 62 * mm, 2.2 * mm, stroke=0, fill=1)
        canvas.setFillColor(HexColor("#B8CBDA"))
        canvas.setFont("Sans", 8.5)
        canvas.drawString(18 * mm, 19 * mm, "Deep-research report | Evidence checked against primary sources")
    else:
        canvas.setStrokeColor(RULE)
        canvas.setLineWidth(0.45)
        canvas.line(18 * mm, height - 12 * mm, width - 18 * mm, height - 12 * mm)
        canvas.setFont("Sans-Bold", 7.4)
        canvas.setFillColor(BLUE)
        canvas.drawString(18 * mm, height - 9.2 * mm, "OSCILLONS - RESEARCH ORIENTATION")
        canvas.setFont("Sans", 7.4)
        canvas.setFillColor(MUTED)
        canvas.drawRightString(width - 18 * mm, height - 9.2 * mm, "6 SEPTEMBER 2026")
        canvas.line(18 * mm, 12 * mm, width - 18 * mm, 12 * mm)
        canvas.setFont("Sans", 7.3)
        canvas.drawString(18 * mm, 8.2 * mm, "Big picture | source audit | quantitative roadmap")
        canvas.drawRightString(width - 18 * mm, 8.2 * mm, str(page - 1))
    canvas.restoreState()


def arrow(d: Drawing, x1: float, y1: float, x2: float, y2: float, color=BLUE, width=1.3):
    d.add(Line(x1, y1, x2, y2, strokeColor=color, strokeWidth=width))
    import math

    angle = math.atan2(y2 - y1, x2 - x1)
    size = 5
    p1 = (x2, y2)
    p2 = (x2 - size * math.cos(angle - 0.55), y2 - size * math.sin(angle - 0.55))
    p3 = (x2 - size * math.cos(angle + 0.55), y2 - size * math.sin(angle + 0.55))
    d.add(Polygon([p1[0], p1[1], p2[0], p2[1], p3[0], p3[1]], fillColor=color, strokeColor=color))


def balance_diagram(styles):
    d = Drawing(480, 128)
    cards = [
        (0, PALE_BLUE, BLUE, "DISPERSION", ["Gradients spread", "localized energy"]),
        (177, PALE_TEAL, TEAL, "OSCILLON CORE", ["Attractive nonlinearity", "supports coherence"]),
        (354, PALE_CORAL, CORAL, "RADIATIVE LEAK", ["Harmonics above m", "carry energy away"]),
    ]
    for x, fill, edge, title, lines in cards:
        d.add(Rect(x, 32, 126, 74, rx=8, ry=8, fillColor=fill, strokeColor=edge, strokeWidth=1.1))
        d.add(String(x + 63, 83, title, fontName="Sans-Bold", fontSize=8.3, textAnchor="middle", fillColor=edge))
        d.add(String(x + 63, 63, lines[0], fontName="Sans", fontSize=7.4, textAnchor="middle", fillColor=INK))
        d.add(String(x + 63, 51, lines[1], fontName="Sans", fontSize=7.4, textAnchor="middle", fillColor=INK))
    arrow(d, 128, 69, 173, 69, BLUE)
    arrow(d, 303, 69, 350, 69, CORAL)
    d.add(String(151, 78, "balance", fontName="Sans-Italic", fontSize=6.8, textAnchor="middle", fillColor=MUTED))
    d.add(String(327, 78, "slow", fontName="Sans-Italic", fontSize=6.8, textAnchor="middle", fillColor=MUTED))
    caption = Paragraph(
        "Figure 1. Localization is a dynamical balance, not an exact prohibition on decay.", styles["Caption"]
    )
    return [Spacer(1, 4), d, caption]


def hierarchy_diagram(styles):
    d = Drawing(480, 206)
    boxes = [
        (120, 158, 240, 34, NAVY, colors.white, "FULL REAL KLEIN-GORDON FIELD", "radiation, relativistic harmonics, exact energy"),
        (15, 94, 210, 40, BLUE, colors.white, "QB / PQB / HARMONIC DESCRIPTION", "instantaneous core, outgoing-tail matching"),
        (255, 94, 210, 40, TEAL, colors.white, "NONRELATIVISTIC EFFECTIVE FIELD", "gradient expansion, approximate U(1)"),
        (255, 28, 210, 38, GOLD, NAVY, "CUBIC NLS TRUNCATION", "useful scaling model; narrowest domain"),
    ]
    for x, y, w, h, fill, fg, title, subtitle in boxes:
        d.add(Rect(x, y, w, h, rx=7, ry=7, fillColor=fill, strokeColor=fill))
        d.add(String(x + w / 2, y + h - 13, title, fontName="Sans-Bold", fontSize=7.7, textAnchor="middle", fillColor=fg))
        d.add(String(x + w / 2, y + 9, subtitle, fontName="Sans", fontSize=6.7, textAnchor="middle", fillColor=fg))
    arrow(d, 240, 158, 124, 136, BLUE)
    arrow(d, 240, 158, 360, 136, TEAL)
    arrow(d, 360, 93, 360, 69, GOLD)
    d.add(String(76, 145, "Fourier / boundary construction", fontName="Sans-Italic", fontSize=6.5, fillColor=MUTED))
    d.add(String(343, 145, "slow / large regime", fontName="Sans-Italic", fontSize=6.5, fillColor=MUTED))
    d.add(String(373, 78, "truncate", fontName="Sans-Italic", fontSize=6.5, fillColor=MUTED))
    d.add(Rect(15, 3, 450, 16, rx=4, ry=4, fillColor=PALE_CORAL, strokeColor=CORAL, strokeWidth=0.7))
    d.add(String(240, 8, "A singularity in a lower box is not automatically a singularity in a higher box.", fontName="Sans-Bold", fontSize=7.1, textAnchor="middle", fillColor=CORAL))
    caption = Paragraph("Figure 2. Model hierarchy and the direction in which assumptions accumulate.", styles["Caption"])
    return [Spacer(1, 4), d, caption]


def validation_diagram(styles):
    d = Drawing(480, 184)
    labels = [
        ("1", "LINEAR BASELINE", "exact or manufactured solution", PALE_BLUE, BLUE),
        ("2", "ORDER TEST", "refine space and time separately", PALE_TEAL, TEAL),
        ("3", "BOUNDARY TEST", "move wall and vary absorber", PALE_GOLD, HexColor("#B28218")),
        ("4", "BOUNDED BENCHMARK", "reproduce one literature trajectory", PALE_BLUE, BLUE),
        ("5", "PHYSICS SCAN", "only now classify formation and decay", PALE_CORAL, CORAL),
    ]
    y = 145
    for idx, title, sub, fill, edge in labels:
        d.add(Circle(28, y + 14, 12, fillColor=edge, strokeColor=edge))
        d.add(String(28, y + 10.5, idx, fontName="Sans-Bold", fontSize=9, textAnchor="middle", fillColor=colors.white))
        d.add(Rect(52, y, 398, 29, rx=5, ry=5, fillColor=fill, strokeColor=edge, strokeWidth=0.8))
        d.add(String(67, y + 17, title, fontName="Sans-Bold", fontSize=7.9, fillColor=edge))
        d.add(String(195, y + 17, sub, fontName="Sans", fontSize=7.4, fillColor=INK))
        if idx != "5":
            arrow(d, 28, y - 2, 28, y - 14, MUTED, 0.9)
        y -= 35
    caption = Paragraph("Figure 3. Each rung removes a different false-positive mechanism.", styles["Caption"])
    return [Spacer(1, 4), d, caption]


def callout(text: str, styles) -> Table:
    p = paragraph(text, styles["Callout"])
    table = Table([[p]], colWidths=[174 * mm], hAlign="LEFT")
    table.setStyle(
        TableStyle(
            [
                ("BACKGROUND", (0, 0), (-1, -1), PALE_GOLD),
                ("BOX", (0, 0), (-1, -1), 0.5, HexColor("#E3D29C")),
                ("LINEBEFORE", (0, 0), (0, -1), 4, GOLD),
                ("LEFTPADDING", (0, 0), (-1, -1), 11),
                ("RIGHTPADDING", (0, 0), (-1, -1), 10),
                ("TOPPADDING", (0, 0), (-1, -1), 8),
                ("BOTTOMPADDING", (0, 0), (-1, -1), 8),
            ]
        )
    )
    table.spaceBefore = 4
    table.spaceAfter = 9
    return table


def equation_box(lines: list[str], styles) -> Table:
    content = "<br/>".join(html.escape(line.strip(), quote=False) for line in lines if line.strip())
    p = Paragraph(content, styles["Equation"])
    t = Table([[p]], colWidths=[164 * mm], hAlign="CENTER")
    t.setStyle(
        TableStyle(
            [
                ("BACKGROUND", (0, 0), (-1, -1), PALE_BLUE),
                ("BOX", (0, 0), (-1, -1), 0.55, RULE),
                ("LEFTPADDING", (0, 0), (-1, -1), 9),
                ("RIGHTPADDING", (0, 0), (-1, -1), 9),
                ("TOPPADDING", (0, 0), (-1, -1), 6),
                ("BOTTOMPADDING", (0, 0), (-1, -1), 6),
            ]
        )
    )
    t.spaceBefore = 3
    t.spaceAfter = 8
    return t


def markdown_table(rows: list[list[str]], styles) -> Table:
    ncols = len(rows[0])
    char_weights = []
    for col in range(ncols):
        lengths = [len(re.sub(r"\[[^\]]+\]\([^)]+\)", "link", row[col])) for row in rows]
        char_weights.append(max(10, min(43, max(lengths))))
    total_width = 174 * mm
    weight_sum = sum(char_weights)
    widths = [total_width * w / weight_sum for w in char_weights]
    min_width = 24 * mm if ncols <= 4 else 18 * mm
    widths = [max(min_width, w) for w in widths]
    widths = [w * total_width / sum(widths) for w in widths]

    data = []
    for r_idx, row in enumerate(rows):
        style = styles["TableHead"] if r_idx == 0 else styles["Table"]
        data.append([paragraph(cell.strip(), style) for cell in row])
    t = Table(data, colWidths=widths, repeatRows=1, hAlign="LEFT", splitByRow=1)
    commands = [
        ("BACKGROUND", (0, 0), (-1, 0), NAVY),
        ("TEXTCOLOR", (0, 0), (-1, 0), colors.white),
        ("GRID", (0, 0), (-1, -1), 0.35, RULE),
        ("VALIGN", (0, 0), (-1, -1), "TOP"),
        ("LEFTPADDING", (0, 0), (-1, -1), 5),
        ("RIGHTPADDING", (0, 0), (-1, -1), 5),
        ("TOPPADDING", (0, 0), (-1, -1), 4.5),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 4.5),
    ]
    for r_idx in range(1, len(data)):
        if r_idx % 2 == 0:
            commands.append(("BACKGROUND", (0, r_idx), (-1, r_idx), HexColor("#F6F8FA")))
    t.setStyle(TableStyle(commands))
    t.spaceBefore = 4
    t.spaceAfter = 8
    return t


def parse_markdown(md: str, styles):
    lines = md.splitlines()
    body_start = lines.index("<!-- REPORT_BODY -->") + 1
    lines = lines[body_start:]
    story = []

    # Cover page
    story.extend(
        [
            Spacer(1, 66 * mm),
            Paragraph("DEEP-RESEARCH ORIENTATION", styles["CoverKicker"]),
            Paragraph("Oscillons", styles["CoverTitle"]),
            Paragraph(
                "Big picture, trustworthy foundations, and a quantitative path through the material in this folder",
                styles["CoverSub"],
            ),
            HRFlowable(width="36%", thickness=2, color=TEAL, hAlign="LEFT", spaceBefore=2, spaceAfter=13),
            Paragraph(
                "Prepared for an early-stage research project<br/>Real scalar fields | radiation | dimensional dependence | terminal decay | radial numerics",
                styles["CoverMeta"],
            ),
            PageBreak(),
            Paragraph("Executive orientation", styles["Executive"]),
        ]
    )

    i = 0
    while i < len(lines):
        raw = lines[i]
        line = raw.strip()
        if not line:
            i += 1
            continue
        if line == "<!-- PAGEBREAK -->":
            story.append(PageBreak())
            i += 1
            continue
        if line == "[[TOC]]":
            story.append(Spacer(1, 4))
            story.append(Paragraph("Contents", styles["Heading1"]))
            toc = TableOfContents()
            toc.levelStyles = [
                ParagraphStyle(
                    "TOC1", fontName="Sans-Bold", fontSize=9.4, leading=12, leftIndent=0, firstLineIndent=0, textColor=NAVY, spaceBefore=3
                ),
                ParagraphStyle(
                    "TOC2", fontName="Sans", fontSize=8.2, leading=10.2, leftIndent=14, firstLineIndent=0, textColor=MUTED, spaceBefore=1
                ),
            ]
            story.append(toc)
            i += 1
            continue
        if line == "[[DIAGRAM_BALANCE]]":
            story.extend(balance_diagram(styles))
            i += 1
            continue
        if line == "[[DIAGRAM_HIERARCHY]]":
            story.extend(hierarchy_diagram(styles))
            i += 1
            continue
        if line == "[[DIAGRAM_LADDER]]":
            story.extend(validation_diagram(styles))
            i += 1
            continue
        if line.startswith("### "):
            story.append(paragraph(line[4:], styles["Heading2"]))
            i += 1
            continue
        if line.startswith("## "):
            story.append(paragraph(line[3:], styles["Heading1"]))
            i += 1
            continue
        if line.startswith("# "):
            i += 1
            continue
        if line == "$$":
            eq_lines = []
            i += 1
            while i < len(lines) and lines[i].strip() != "$$":
                eq_lines.append(lines[i])
                i += 1
            story.append(equation_box(eq_lines, styles))
            i += 1
            continue
        if line.startswith("> "):
            quote_lines = []
            while i < len(lines) and lines[i].strip().startswith("> "):
                quote_lines.append(lines[i].strip()[2:])
                i += 1
            story.append(callout(" ".join(quote_lines), styles))
            continue
        if line.startswith("|") and line.endswith("|"):
            table_lines = []
            while i < len(lines):
                candidate = lines[i].strip()
                if not (candidate.startswith("|") and candidate.endswith("|")):
                    break
                table_lines.append(candidate)
                i += 1
            rows = [[c.strip() for c in t.strip("|").split("|")] for t in table_lines]
            if len(rows) >= 2 and all(re.fullmatch(r":?-{3,}:?", c.replace(" ", "")) for c in rows[1]):
                rows.pop(1)
            story.append(markdown_table(rows, styles))
            continue
        if re.match(r"^-\s+", line):
            items = []
            while i < len(lines) and re.match(r"^-\s+", lines[i].strip()):
                item_text = re.sub(r"^-\s+", "", lines[i].strip())
                items.append(ListItem(paragraph(item_text, styles["List"]), leftIndent=12, bulletColor=TEAL))
                i += 1
            story.append(
                ListFlowable(items, bulletType="bullet", start="circle", leftIndent=17, bulletFontName="Sans", bulletFontSize=6.5, spaceAfter=5)
            )
            continue
        if re.match(r"^\d+\.\s+", line):
            items = []
            first_number = int(re.match(r"^(\d+)\.", line).group(1))
            while i < len(lines) and re.match(r"^\d+\.\s+", lines[i].strip()):
                item_text = re.sub(r"^\d+\.\s+", "", lines[i].strip())
                items.append(ListItem(paragraph(item_text, styles["List"]), leftIndent=16))
                i += 1
            story.append(
                ListFlowable(
                    items,
                    bulletType="1",
                    start=str(first_number),
                    leftIndent=20,
                    bulletFontName="Sans-Bold",
                    bulletFontSize=8.5,
                    bulletColor=BLUE,
                    spaceAfter=6,
                )
            )
            continue

        para_lines = [line]
        i += 1
        while i < len(lines):
            nxt = lines[i].strip()
            if not nxt:
                break
            if (
                nxt.startswith("#")
                or nxt.startswith("> ")
                or nxt.startswith("|")
                or nxt == "$$"
                or nxt.startswith("[[")
                or nxt.startswith("<!--")
                or re.match(r"^-\s+", nxt)
                or re.match(r"^\d+\.\s+", nxt)
            ):
                break
            para_lines.append(nxt)
            i += 1
        story.append(paragraph(" ".join(para_lines), styles["Body"]))

    return story


def main() -> None:
    register_fonts()
    styles = make_styles()
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    doc = ResearchDocTemplate(
        str(OUTPUT),
        styles,
        pagesize=A4,
        leftMargin=18 * mm,
        rightMargin=18 * mm,
        topMargin=18 * mm,
        bottomMargin=17 * mm,
        title="Oscillons: A Research Orientation and Simulation Reading Map",
        author="OpenAI Codex",
        subject="Deep-research orientation based on the local oscillon notes and simulation files",
    )
    frame = Frame(doc.leftMargin, doc.bottomMargin, doc.width, doc.height, id="normal")
    doc.addPageTemplates([PageTemplate(id="research", frames=[frame], onPage=cover_background)])
    story = parse_markdown(SOURCE.read_text(encoding="utf-8"), styles)
    doc.multiBuild(story)
    print(OUTPUT)


if __name__ == "__main__":
    main()
