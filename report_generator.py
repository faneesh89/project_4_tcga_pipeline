# report_generator.py
# ─────────────────────────────────────────────────────────────────────────────
# Takes the finished per-gene frequency results and writes an HTML report.
#
# Presentation only — this file doesn't know anything about MAF files,
# mutation biology, or cancer genes. It just formats finished numbers.
# Same structural pattern as Project 3's report_generator.py.
# ─────────────────────────────────────────────────────────────────────────────

import os
from datetime import datetime

from frequency_calculator import format_mutation_types


def _frequency_css_class(freq_pct):
    """Color-code rows by how frequently the gene is mutated."""
    if freq_pct >= 40:
        return 'very-high'
    elif freq_pct >= 20:
        return 'high'
    elif freq_pct >= 10:
        return 'moderate'
    else:
        return 'low'


def generate_html_report(results, stats, output_file="output/cancer_report.html",
                          total_patients=None):
    """
    results: list of per-gene dicts from frequency_calculator.calculate_frequencies()
    stats:   the stats dict from mutation_analyzer.analyze_mutations()
    output_file: where to write the HTML
    """

    # use the TRUE total passed from main.py; fall back only if absent
    if total_patients is None:
        total_patients = len(stats.get('all_patients', set()))
    patients_with_drivers = len(stats.get('all_patients', set()))
    total_mutations = stats.get('total_mutations', 0)
    severity = stats.get('severity_counts', {})

    rows_html = ""
    for r in results:
        css_class = _frequency_css_class(r['frequency_pct'])
        types_str = format_mutation_types(r['mutation_types'])
        rows_html += f"""
        <tr class="{css_class}">
            <td><strong>{r['gene']}</strong></td>
            <td>{r['patient_count']} / {r['total_patients']}</td>
            <td><strong>{r['frequency_pct']}%</strong></td>
            <td>{r['mutation_count']}</td>
            <td>{types_str}</td>
        </tr>
        """

    if not rows_html:
        rows_html = ('<tr><td colspan="5">No driver gene mutations '
                     'found in this dataset.</td></tr>')

    severity_html = ""
    for level in ['high', 'moderate', 'low', 'other']:
        count = severity.get(level, 0)
        if count:
            severity_html += f'<span class="sev sev-{level}">{level}: {count}</span> '

    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<title>Cancer Driver Gene Mutation Report</title>
<style>
    body {{
        font-family: 'Segoe UI', Arial, sans-serif;
        background: #f4f6f8;
        color: #2c3e50;
        margin: 0;
        padding: 30px;
    }}
    .container {{
        max-width: 1000px;
        margin: 0 auto;
        background: white;
        border-radius: 8px;
        box-shadow: 0 2px 10px rgba(0,0,0,0.08);
        overflow: hidden;
    }}
    .header {{
        background: #7B1FA2;
        color: white;
        padding: 25px 30px;
    }}
    .header h1 {{ margin: 0 0 5px 0; font-size: 24px; }}
    .header p {{ margin: 0; color: #E1BEE7; font-size: 13px; }}
    .content {{ padding: 25px 30px; }}
    .summary {{
        background: #F3E5F5;
        border-left: 4px solid #7B1FA2;
        padding: 15px;
        margin-bottom: 20px;
        font-size: 14px;
    }}
    table {{ width: 100%; border-collapse: collapse; margin-top: 10px; }}
    th {{
        background: #9C27B0;
        color: white;
        text-align: left;
        padding: 10px 12px;
        font-size: 13px;
    }}
    td {{
        padding: 10px 12px;
        border-bottom: 1px solid #e0e0e0;
        font-size: 13px;
        vertical-align: top;
    }}
    tr.very-high {{ background: #FFEBEE; }}
    tr.high {{ background: #FFF3E0; }}
    tr.moderate {{ background: #FFFDE7; }}
    tr.low {{ background: #F5F5F5; }}
    .sev {{
        display: inline-block;
        padding: 4px 12px;
        border-radius: 3px;
        margin-right: 8px;
        font-size: 12px;
    }}
    .sev-high {{ background: #FFCDD2; }}
    .sev-moderate {{ background: #FFE0B2; }}
    .sev-low {{ background: #DCEDC8; }}
    .sev-other {{ background: #E0E0E0; }}
    .disclaimer {{
        margin-top: 25px;
        padding: 15px;
        background: #FFF3E0;
        border-left: 4px solid #E65100;
        font-size: 12px;
        color: #555;
    }}
    .footer {{
        text-align: center;
        padding: 15px;
        font-size: 11px;
        color: #999;
        border-top: 1px solid #eee;
    }}
</style>
</head>
<body>
    <div class="container">
        <div class="header">
            <h1>Cancer Driver Gene Mutation Report</h1>
            <p>TCGA somatic mutation analysis &middot;
               Generated {datetime.now().strftime('%B %d, %Y at %H:%M')}</p>
        </div>
        <div class="content">
            <div class="summary">
                <strong>{total_patients}</strong> patients analyzed &middot;
                <strong>{patients_with_drivers}</strong> with driver mutations &middot;
                <strong>{total_mutations}</strong> driver-gene mutations &middot;
                <strong>{len(results)}</strong> genes affected
                <br><br>
                Mutation severity breakdown: {severity_html or '—'}
            </div>

            <table>
                <tr>
                    <th>Gene</th>
                    <th>Patients Mutated</th>
                    <th>Frequency</th>
                    <th>Total Mutations</th>
                    <th>Mutation Types</th>
                </tr>
                {rows_html}
            </table>

            <div class="disclaimer">
                <strong>Note:</strong> This analysis uses a curated subset of
                well-documented cancer driver genes, not the complete COSMIC
                Cancer Gene Census (719 genes). "Patients Mutated" counts unique
                patients per gene &mdash; a patient with multiple mutations in the
                same gene is counted once. Frequency percentages reflect this
                cohort only and should not be generalized. This is an
                educational/portfolio project, not a clinical diagnostic tool.
            </div>
        </div>
        <div class="footer">
            TCGA Cancer Variant Analysis Pipeline &mdash; Project 4
        </div>
    </div>
</body>
</html>"""

    out_dir = os.path.dirname(output_file)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    with open(output_file, 'w', encoding='utf-8') as f:
        f.write(html)

    print(f"Report saved to {output_file}")
    return output_file