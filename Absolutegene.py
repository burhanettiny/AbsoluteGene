import streamlit as st
import pandas as pd
import numpy as np
import plotly.graph_objects as go
import scipy.stats as stats
from io import BytesIO
from reportlab.lib.pagesizes import letter
from reportlab.lib import colors
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.platypus import (Table, TableStyle, SimpleDocTemplate, Paragraph,
                                 Spacer, PageBreak, Image as RLImage, HRFlowable, Flowable)
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.pdfbase import pdfmetrics
import matplotlib.pyplot as plt
import os

try:
    plt.rcParams['font.family'] = 'DejaVu Sans'
    plt.rcParams['axes.unicode_minus'] = False
except Exception:
    pass

st.set_page_config(
    page_title="AbsoluteGene",
    page_icon="🧪",
    layout="wide"
)

# ═══════════════════════════════════════════════════════════════════════════════
# GLOBAL CSS
# ═══════════════════════════════════════════════════════════════════════════════
st.markdown("""
<style>
.block-container { padding-top: 1.2rem !important; padding-bottom: 1rem !important; }
div[data-testid="stAlert"] { padding: 5px 10px !important; font-size: 12px !important; }
div[data-testid="stNumberInput"] { margin-bottom: 0 !important; }
div[data-testid="stRadio"]       { margin-bottom: 0 !important; }
section[data-testid="stSidebar"] > div:first-child { padding-top: 0.8rem !important; }
section[data-testid="stSidebar"] hr { margin: 4px 0 !important; }
</style>
""", unsafe_allow_html=True)

if 'language' not in st.session_state:
    st.session_state.language = "English"

flags = {"Türkçe": "🇹🇷", "English": "🇬🇧"}
default_index = list(flags.keys()).index(st.session_state.language) if st.session_state.language in flags else 1

selected_language = st.sidebar.selectbox(
    "🌐 Language / Dil",
    options=[f"{flags[lang]} {lang}" for lang in flags],
    index=default_index,
    label_visibility="collapsed"
)
selected_language_name = selected_language.split(' ', 1)[1]
language_map = {"Türkçe": "tr", "English": "en"}
language_code = language_map.get(selected_language_name, "en")

# ═══════════════════════════════════════════════════════════════════════════════
# TRANSLATIONS (TR / EN)
# ═══════════════════════════════════════════════════════════════════════════════
translations = {
    "tr": {
        "title": "🧪 AbsoluteGene: Dijital PCR (dPCR/ddPCR) Gen Ekspresyonu ve CNV Analizi",
        "subtitle": "B. Yalçınkaya tarafından geliştirildi — GeneQuantify'nin dijital PCR versiyonu",
        "tab_data": "Veri Girişi",
        "tab_results": "Sonuçlar",
        "tab_report": "Rapor",
        "study_design": "⚙️ Çalışma Tasarımı",
        "num_target_genes": "🔹 Hedef Gen Sayısı",
        "num_patient_groups": "🔹 Hasta Grubu Sayısı",
        "num_ref_genes": "🔹 Referans Gen Sayısı",
        "ref_gene_help": "dMIQE kılavuzları sağlam normalizasyon için ≥1 doğrulanmış referans lokusu önerir; ≥2 önerilir.",
        "ploidy_label": "🔹 Referans Lokus Ploidi Sayısı",
        "ploidy_help": "Kopya sayısı varyasyonu (CNV) hesabı için referans lokusun bilinen kopya sayısı (diploid organizmada genelde 2).",
        "partition_vol_label": "🔹 Partisyon Hacmi (nL)",
        "partition_vol_help": "Damlacık/kuyucuk başına hacim. Bio-Rad QX200: ~0.85 nL, QX ONE: ~0.7 nL, Qiagen QIAcuity: plaka tipine göre değişir. Sadece kopya/µL hesaplaması için kullanılır; oran/CNV hesabını etkilemez.",
        "qc_min_partitions": "🔹 Minimum Kabul Edilebilir Partisyon Sayısı (QC)",
        "qc_min_partitions_help": "Bu değerin altındaki replikatlar düşük kalite olarak işaretlenir (dMIQE önerisi: ddPCR için tipik olarak ≥10.000).",
        "outlier_section_title": "### 🔍 Aykırı Değer Tespiti Ayarları",
        "outlier_enable": "Aykırı değer tespitini etkinleştir",
        "outlier_enable_help": "λ (partisyon başına kopya sayısı) değerlerindeki istatistiksel olarak aşırı replikatları tespit eder.",
        "outlier_method_label": "Tespit yöntemi",
        "outlier_method_help": "Grubbs: normal dağılım için. IQR: parametrik olmayan, çarpık dağılımlar için.",
        "outlier_alpha_label": "Anlamlılık düzeyi (α)",
        "outlier_iqr_label": "IQR çarpanı (k)",
        "patient_data_header": "📥 Hasta ve Kontrol Grubu Partisyon Verisi Girin",
        "target_gene": "Hedef Gen",
        "reference_gene": "Referans Gen",
        "control_group": "🧬 Kontrol Grubu",
        "patient_group": "🩸 Hasta Grubu",
        "positive_partitions": "Pozitif Partisyon Sayısı",
        "total_partitions": "Toplam Kabul Edilen Partisyon Sayısı",
        "input_format_info": "ℹ️ Her satıra bir replikat gelecek şekilde girin. 'Pozitif' ve 'Toplam' kutularındaki satır sayıları eşleşmelidir.",
        "warning_empty_input": "⚠️ Dikkat: Verileri alt alta yazın; Pozitif ve Toplam kutularındaki satır sayısı eşit olmalı.",
        "qc_fail_warning": "⚠️ **Düşük partisyon sayısı uyarısı:** {n} replikat, minimum eşik ({thr}) altında kabul edilen partisyona sahip. Bu replikatlar analiz için düşük güvenilirlikte olabilir.",
        "saturation_warning": "⚠️ **Doygunluk uyarısı:** {n} replikatta tüm partisyonlar pozitif (p≥1.0) — λ hesaplanamıyor. Örneğin daha fazla seyreltilmesi gerekiyor.",
        "gr_tbl": "📋 Giriş Verileri Tablosu (λ ve %95 GA dahil)",
        "nil_mine": "📊 Sonuçlar",
        "download_csv": "📥 CSV İndir",
        "generate_pdf": "📥 PDF Raporu Hazırla",
        "pdf_report": "Dijital PCR Analiz Raporu",
        "sample_number": "Replikat No",
        "lambda_col": "λ (kopya/partisyon)",
        "ci_low_col": "%95 GA Alt",
        "ci_high_col": "%95 GA Üst",
        "conc_col": "Konsantrasyon (kopya/µL)",
        "qc_col": "QC Durumu",
        "qc_pass": "✅ Geçti",
        "qc_fail": "❌ Düşük n",
        "qc_saturated": "❌ Doygun",
        "ratio_col": "Normalize Oran (Hedef/Referans)",
        "cn_col": "Kopya Sayısı (CN)",
        "fc_col": "Kat Değişimi (vs Kontrol)",
        "regulation_status": "Regülasyon Durumu",
        "no_change": "Değişim Yok",
        "upregulated": "Yukarı Regüle / Kazanım",
        "downregulated": "Aşağı Regüle / Kayıp",
        "outlier_excluded_no": "Hayır",
        "outlier_excluded_yes": "Evet",
        "genorm_title": "Referans Lokus Stabilitesi",
        "stable": "Kararlı",
        "borderline": "Sınırda",
        "unstable": "Kararsız",
        "m_value": "M-değeri",
        "method_comparison": "📊 Sonuç Özeti",
        "statistical_results": "📈 İstatistiksel Sonuçlar",
        "test_type": "Test Türü",
        "test_method": "Kullanılan Test",
        "test_pvalue": "Test P-değeri",
        "significance": "Anlamlılık",
        "significant": "Anlamlı",
        "insignificant": "Anlamsız",
        "parametric": "Parametrik",
        "non_parametric": "Nonparametrik",
        "t_test": "t-test",
        "welch_t_test": "Welch t-testi",
        "mann_whitney_u_test": "Mann-Whitney U testi",
        "stat_decision_title": "🔬 İstatistiksel karar",
        "stat_decision_steps": "**Adım adım test seçimi:**",
        "stat_shapiro_title": "**1. Shapiro-Wilk normallik testi**",
        "stat_normal": "Normal",
        "stat_nonnormal": "Normal değil",
        "stat_levene_title": "**2. Levene varyans homojenliği testi**",
        "stat_levene_skipped": "**2. Levene testi** — *atlandı* (normallik sağlanmadı)",
        "stat_equal_var": "Eşit varyans",
        "stat_unequal_var": "Eşitsiz varyans",
        "stat_selected_test": "**3. Seçilen test:**",
        "stat_reason": "**Gerekçe:**",
        "stat_result": "**Sonuç:**",
        "stat_reason_nonnormal": "Bir veya her iki grupta normal dağılım sağlanmadı",
        "stat_reason_normal_equal": "Her iki grup normal + eşit varyans",
        "stat_reason_normal_unequal": "Her iki grup normal + eşitsiz varyans",
        "stat_multigroup_note": "⚠️ Not: ≥3 grup varsa, aşağıdaki Çoklu Grup Karşılaştırması bölümüne bakın.",
        "multigroup_title": "## 📊 Çoklu Grup Karşılaştırma Analizi",
        "multigroup_omnibus_test": "Omnibus Testi",
        "multigroup_pvalue": "p-değeri",
        "multigroup_result": "Sonuç",
        "multigroup_significant": "Anlamlı",
        "multigroup_not_significant": "Anlamlı değil",
        "multigroup_omnibus_ns": "ℹ️ Omnibus testi anlamlı değil (p ≥ 0.05). Post-hoc karşılaştırmalar bilgi amaçlıdır.",
        "multigroup_posthoc_label": "**Post-hoc:**",
        "multigroup_dl_button": "📥 Post-hoc sonuçlarını indir —",
        "multigroup_2group_note": "ℹ️ Yalnızca 2 grup tespit edildi. İkili istatistikler yukarıda raporlanmıştır.",
        "multigroup_decision_normal_equal": "✅ Normal + eşit varyans → **Tek yönlü ANOVA + Tukey HSD**",
        "multigroup_decision_normal_unequal": "⚠️ Normal + eşitsiz varyans → **Welch ANOVA + Games-Howell**",
        "multigroup_decision_nonnormal": "⚠️ Normal değil → **Kruskal-Wallis + Dunn post-hoc**",
        "multigene_title": "### 🧬 Çoklu Gen Çoklu Karşılaştırma Düzeltmesi",
        "multigene_sig_raw": "Anlamlı (ham)",
        "multigene_sig_bonf": "Anlamlı (Bonferroni)",
        "multigene_sig_fdr": "Anlamlı (FDR B-H)",
        "multigene_warning": "⚠️ Düzeltme sonrası {lost} sonuç FDR düzeltmesi sonrası anlamlılığını yitirdi.",
        "multigene_success": "✅ {n} anlamlı sonucun tamamı FDR düzeltmesi sonrasında da anlamlı kalmaktadır.",
        "multigene_no_sig": "Ham p < 0.05 eşiğinde anlamlı ikili sonuç tespit edilmedi.",
        "multigene_dl_button": "📥 Düzeltilmiş p-değerlerini indir (CSV)",
        "multigene_1gene_note": "ℹ️ Yalnızca 1 hedef gen analiz edildi — genler arası düzeltme uygulanamaz.",
        "distribution_graph": "Dağılım Grafiği",
        "x_axis_title": "Grup Adı",
        "dist_plot_mode_label": "📊 Dağılım Grafiği — Görüntüleme Modu",
        "dist_plot_ratio": "Normalize Oran (Hedef/Referans) — önerilen",
        "dist_plot_lambda": "λ (kopya/partisyon) — ham",
        "dist_plot_fc": "Kat Değişimi — kontrole göre",
        "error_no_data": "Veri bulunamadı, PDF oluşturulamadı.",
        "pdf_ready": "{n} kayıt hazır — PDF oluşturabilirsiniz.",
        "sidebar_github_btn": "⭐ GitHub'da Kaynak Kodu Görüntüle",
        "sidebar_sister_tool": "🧬 qPCR için: GeneQuantify",
        "rdml_expander": "ℹ️ Veri girişi hakkında",
        "guide_btn": "📘 Kullanım Kılavuzu",
    },
    "en": {
        "title": "🧪 AbsoluteGene: Digital PCR (dPCR/ddPCR) Gene Expression & CNV Analysis",
        "subtitle": "Developed by B. Yalçınkaya — the digital PCR companion to GeneQuantify",
        "tab_data": "Data Entry",
        "tab_results": "Results",
        "tab_report": "Report",
        "study_design": "⚙️ Study Design",
        "num_target_genes": "🔹 Number of Target Genes",
        "num_patient_groups": "🔹 Number of Patient Groups",
        "num_ref_genes": "🔹 Number of Reference Genes",
        "ref_gene_help": "dMIQE guidelines recommend ≥1 validated reference locus for normalization; ≥2 is preferred.",
        "ploidy_label": "🔹 Reference Locus Ploidy",
        "ploidy_help": "Known copy number of the reference locus, used for copy number variation (CNV) calculation (typically 2 for a diploid organism).",
        "partition_vol_label": "🔹 Partition Volume (nL)",
        "partition_vol_help": "Volume per droplet/well. Bio-Rad QX200: ~0.85 nL, QX ONE: ~0.7 nL, Qiagen QIAcuity: varies by plate type. Only affects copies/µL, not the ratio/CNV calculation.",
        "qc_min_partitions": "🔹 Minimum Accepted Partitions (QC)",
        "qc_min_partitions_help": "Replicates below this value are flagged as low quality (dMIQE recommendation: typically ≥10,000 for ddPCR).",
        "outlier_section_title": "### 🔍 Outlier Detection Settings",
        "outlier_enable": "Enable outlier detection",
        "outlier_enable_help": "Detects statistically extreme replicates in λ (copies per partition).",
        "outlier_method_label": "Detection method",
        "outlier_method_help": "Grubbs: for normally distributed data. IQR: non-parametric, robust for skewed distributions.",
        "outlier_alpha_label": "Significance level (α)",
        "outlier_iqr_label": "IQR multiplier (k)",
        "patient_data_header": "📥 Enter Patient and Control Group Partition Data",
        "target_gene": "Target Gene",
        "reference_gene": "Reference Gene",
        "control_group": "🧬 Control Group",
        "patient_group": "🩸 Patient Group",
        "positive_partitions": "Positive Partition Count",
        "total_partitions": "Total Accepted Partition Count",
        "input_format_info": "ℹ️ Enter one replicate per line. The number of lines in 'Positive' and 'Total' boxes must match.",
        "warning_empty_input": "⚠️ Warning: Enter data one value per line; the number of lines in Positive and Total boxes must be equal.",
        "qc_fail_warning": "⚠️ **Low partition count warning:** {n} replicate(s) have accepted partitions below the minimum threshold ({thr}). These replicates may have reduced reliability.",
        "saturation_warning": "⚠️ **Saturation warning:** {n} replicate(s) have all partitions positive (p≥1.0) — λ cannot be calculated. Further dilution is recommended.",
        "gr_tbl": "📋 Input Data Table (incl. λ and 95% CI)",
        "nil_mine": "📊 Results",
        "download_csv": "📥 Download CSV",
        "generate_pdf": "📥 Prepare PDF Report",
        "pdf_report": "Digital PCR Analysis Report",
        "sample_number": "Replicate #",
        "lambda_col": "λ (copies/partition)",
        "ci_low_col": "95% CI Lower",
        "ci_high_col": "95% CI Upper",
        "conc_col": "Concentration (copies/µL)",
        "qc_col": "QC Status",
        "qc_pass": "✅ Passed",
        "qc_fail": "❌ Low n",
        "qc_saturated": "❌ Saturated",
        "ratio_col": "Normalized Ratio (Target/Reference)",
        "cn_col": "Copy Number (CN)",
        "fc_col": "Fold Change (vs Control)",
        "regulation_status": "Regulation Status",
        "no_change": "No Change",
        "upregulated": "Upregulated / Gain",
        "downregulated": "Downregulated / Loss",
        "outlier_excluded_no": "No",
        "outlier_excluded_yes": "Yes",
        "genorm_title": "Reference Locus Stability",
        "stable": "Stable",
        "borderline": "Borderline",
        "unstable": "Unstable",
        "m_value": "M-value",
        "method_comparison": "📊 Result Summary",
        "statistical_results": "📈 Statistical Results",
        "test_type": "Test Type",
        "test_method": "Test Method",
        "test_pvalue": "Test P-value",
        "significance": "Significance",
        "significant": "Significant",
        "insignificant": "Insignificant",
        "parametric": "Parametric",
        "non_parametric": "Nonparametric",
        "t_test": "t-test",
        "welch_t_test": "Welch t-test",
        "mann_whitney_u_test": "Mann-Whitney U test",
        "stat_decision_title": "🔬 Statistical decision",
        "stat_decision_steps": "**Step-by-step test selection:**",
        "stat_shapiro_title": "**1. Shapiro-Wilk normality test**",
        "stat_normal": "Normal",
        "stat_nonnormal": "Non-normal",
        "stat_levene_title": "**2. Levene variance homogeneity test**",
        "stat_levene_skipped": "**2. Levene test** — *skipped* (normality not met)",
        "stat_equal_var": "Equal variances",
        "stat_unequal_var": "Unequal variances",
        "stat_selected_test": "**3. Selected test:**",
        "stat_reason": "**Reason:**",
        "stat_result": "**Result:**",
        "stat_reason_nonnormal": "Non-normal distribution in one or both groups",
        "stat_reason_normal_equal": "Both groups normal + equal variances",
        "stat_reason_normal_unequal": "Both groups normal + unequal variances",
        "stat_multigroup_note": "⚠️ Note: When ≥3 groups are present, see the Multi-Group Comparison section below.",
        "multigroup_title": "## 📊 Multi-Group Comparison Analysis",
        "multigroup_omnibus_test": "Omnibus Test",
        "multigroup_pvalue": "p-value",
        "multigroup_result": "Result",
        "multigroup_significant": "Significant",
        "multigroup_not_significant": "Not significant",
        "multigroup_omnibus_ns": "ℹ️ Omnibus test is not significant (p ≥ 0.05). Post-hoc comparisons are shown for completeness.",
        "multigroup_posthoc_label": "**Post-hoc:**",
        "multigroup_dl_button": "📥 Download post-hoc results —",
        "multigroup_2group_note": "ℹ️ Only 2 groups detected. Pairwise statistics are reported above.",
        "multigroup_decision_normal_equal": "✅ Normal + equal variances → **One-way ANOVA + Tukey HSD**",
        "multigroup_decision_normal_unequal": "⚠️ Normal + unequal variances → **Welch ANOVA + Games-Howell**",
        "multigroup_decision_nonnormal": "⚠️ Non-normal → **Kruskal-Wallis + Dunn post-hoc**",
        "multigene_title": "### 🧬 Multi-Gene Multiple Comparison Correction",
        "multigene_sig_raw": "Significant (raw)",
        "multigene_sig_bonf": "Significant (Bonferroni)",
        "multigene_sig_fdr": "Significant (FDR B-H)",
        "multigene_warning": "⚠️ After correction, {lost} result(s) no longer significant after FDR adjustment.",
        "multigene_success": "✅ All {n} significant result(s) remain significant after FDR correction.",
        "multigene_no_sig": "No significant pairwise results detected (raw p < 0.05).",
        "multigene_dl_button": "📥 Download corrected p-values (CSV)",
        "multigene_1gene_note": "ℹ️ Only 1 target gene analysed — multi-gene correction not applicable.",
        "distribution_graph": "Distribution Graph",
        "x_axis_title": "Group Name",
        "dist_plot_mode_label": "📊 Distribution Plot — Display Mode",
        "dist_plot_ratio": "Normalized Ratio (Target/Reference) — recommended",
        "dist_plot_lambda": "λ (copies/partition) — raw",
        "dist_plot_fc": "Fold Change — relative to control",
        "error_no_data": "No data found, PDF could not be generated.",
        "pdf_ready": "{n} records ready — you can generate the PDF.",
        "sidebar_github_btn": "⭐ View Source on GitHub",
        "sidebar_sister_tool": "🧬 For qPCR: GeneQuantify",
        "rdml_expander": "ℹ️ About data entry",
        "guide_btn": "📘 User Guide",
    }
}
_t = translations[language_code]

# ═══════════════════════════════════════════════════════════════════════════════
# PDF FONT SYSTEM
# ═══════════════════════════════════════════════════════════════════════════════
def _find_font(candidates):
    import glob as _glob
    for p in candidates:
        if os.path.exists(p):
            return p
    all_ttf = _glob.glob('/usr/share/fonts/**/*.ttf', recursive=True)
    return all_ttf[0] if all_ttf else None

_NOTO_REGULAR = _find_font([
    '/usr/share/fonts/truetype/noto/NotoSans-Regular.ttf',
    '/usr/share/fonts/truetype/freefont/FreeSans.ttf',
    '/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf',
])
_NOTO_BOLD = _find_font([
    '/usr/share/fonts/truetype/noto/NotoSans-Bold.ttf',
    '/usr/share/fonts/truetype/freefont/FreeSansBold.ttf',
    '/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf',
])
try:
    pdfmetrics.registerFont(TTFont('NotoSans', _NOTO_REGULAR))
    pdfmetrics.registerFont(TTFont('NotoSans-Bold', _NOTO_BOLD))
    PDF_FONT, PDF_FONT_BOLD = 'NotoSans', 'NotoSans-Bold'
except Exception:
    PDF_FONT, PDF_FONT_BOLD = 'Helvetica', 'Helvetica-Bold'

def safe_str(text):
    if not isinstance(text, str):
        text = str(text)
    return text.replace('&', '&amp;').replace('<', '&lt;').replace('>', '&gt;')

# ═══════════════════════════════════════════════════════════════════════════════
# CORE dPCR MATH — Poisson statistics
# ═══════════════════════════════════════════════════════════════════════════════
def poisson_lambda(positive, total):
    """
    Returns (lambda, ci_low, ci_high, status) for a single partition (droplet/well) count.
    lambda = copies per partition, derived from the Poisson distribution.
    status: 'ok', 'saturated' (all partitions positive), or 'invalid'.
    95% CI uses the delta-method normal approximation:
        p = positive/total ; lambda = -ln(1-p)
        Var(lambda) ~= p / (n * (1-p))
    """
    if total is None or total <= 0 or positive is None or positive < 0 or positive > total:
        return None, None, None, "invalid"
    p = positive / total
    if p >= 1.0:
        return None, None, None, "saturated"
    lam = -np.log(1 - p)
    se = np.sqrt(p / (total * (1 - p))) if (total * (1 - p)) > 0 else np.nan
    ci_low = max(lam - 1.96 * se, 0.0) if not np.isnan(se) else None
    ci_high = lam + 1.96 * se if not np.isnan(se) else None
    return lam, ci_low, ci_high, "ok"

def geometric_mean(values):
    values = np.array(values, dtype=float)
    values = values[values > 0]
    if len(values) == 0:
        return np.nan
    return np.exp(np.mean(np.log(values)))

def compute_stability_m(ref_lambda_matrix):
    """geNorm-style stability M-value, computed in log-space on lambda values."""
    n_refs, n_samples = ref_lambda_matrix.shape
    if n_refs < 2:
        return np.array([0.0])
    log_mat = np.log(np.clip(ref_lambda_matrix, 1e-9, None))
    m_values = []
    for i in range(n_refs):
        pairwise_sd = []
        for j in range(n_refs):
            if i == j:
                continue
            ratio = log_mat[i] - log_mat[j]
            pairwise_sd.append(np.std(ratio, ddof=1) if len(ratio) > 1 else 0.0)
        m_values.append(np.mean(pairwise_sd))
    return np.array(m_values)

def parse_input_data(input_data):
    values = [x.replace(",", ".").strip() for x in input_data.split() if x.strip()]
    out = []
    for v in values:
        try:
            out.append(float(v))
        except ValueError:
            continue
    return np.array(out)

# ═══════════════════════════════════════════════════════════════════════════════
# OUTLIER DETECTION (Grubbs / IQR) — identical logic to GeneQuantify
# ═══════════════════════════════════════════════════════════════════════════════
def detect_outliers_grubbs(data, alpha=0.05):
    data = np.array(data, dtype=float)
    n = len(data)
    if n < 3:
        return []
    outlier_indices = []
    working = data.copy()
    original_indices = list(range(n))
    while len(working) >= 3:
        mean_w = np.mean(working)
        std_w = np.std(working, ddof=1)
        if std_w == 0:
            break
        g_vals = np.abs(working - mean_w) / std_w
        max_idx = np.argmax(g_vals)
        G = g_vals[max_idx]
        t_crit = stats.t.ppf(1 - alpha / (2 * len(working)), df=len(working) - 2)
        G_crit = ((len(working) - 1) / np.sqrt(len(working))) * \
                 np.sqrt(t_crit**2 / (len(working) - 2 + t_crit**2))
        if G > G_crit:
            outlier_indices.append(original_indices[max_idx])
            original_indices.pop(max_idx)
            working = np.delete(working, max_idx)
        else:
            break
    return outlier_indices

def detect_outliers_iqr(data, multiplier=1.5):
    data = np.array(data, dtype=float)
    q1, q3 = np.percentile(data, [25, 75])
    iqr = q3 - q1
    lower = q1 - multiplier * iqr
    upper = q3 + multiplier * iqr
    return [i for i, v in enumerate(data) if v < lower or v > upper]

def render_outlier_ui(data, label, key_prefix, method, alpha=0.05, k=1.5):
    data = np.array(data, dtype=float)
    detected = detect_outliers_grubbs(data, alpha) if method == "Grubbs" else detect_outliers_iqr(data, k)
    if not detected:
        return data, []
    st.warning(
        f"⚠️ **Potential outlier(s) in {label}** ({method}): "
        f"replicate(s) **{[i+1 for i in detected]}** — λ values: **{[round(data[i], 4) for i in detected]}**"
    )
    excluded = []
    for idx in detected:
        confirm = st.checkbox(
            f"Exclude replicate {idx+1} (λ = {data[idx]:.4f}) from {label}",
            value=False, key=f"{key_prefix}_excl_{idx}"
        )
        if confirm:
            excluded.append(idx)
    if excluded:
        cleaned = np.delete(data, excluded)
        st.info(f"ℹ️ {len(excluded)} replicate(s) excluded from {label}. Remaining n = {len(cleaned)}.")
        return cleaned, excluded
    return data, []

# ═══════════════════════════════════════════════════════════════════════════════
# HEADER
# ═══════════════════════════════════════════════════════════════════════════════
st.markdown(
    f"""
    <div style="background:linear-gradient(90deg,#004d40,#00796b);
                color:white;padding:10px 18px;border-radius:8px;margin-bottom:8px;">
        <div style="font-size:20px;font-weight:800;">{_t['title']}</div>
        <div style="font-size:11px;opacity:0.8;margin-top:2px;">{_t['subtitle']}</div>
    </div>
    """,
    unsafe_allow_html=True
)

# ═══════════════════════════════════════════════════════════════════════════════
# SIDEBAR — User guide, links
# ═══════════════════════════════════════════════════════════════════════════════
guide_clicked = st.sidebar.button(_t['guide_btn'], use_container_width=True)
if guide_clicked:
    @st.dialog("📘 AbsoluteGene — User Guide" if language_code == "en" else "📘 AbsoluteGene — Kullanım Kılavuzu", width="large")
    def show_guide():
        if language_code == "en":
            st.markdown("""
### How dPCR analysis differs from qPCR (GeneQuantify)

Digital PCR partitions each sample into thousands of micro-reactions (droplets or
wells). Each partition ends up **positive** (target present) or **negative**
(target absent). Instead of a Cq value, the input here is:

- **Positive partition count** (per replicate)
- **Total accepted partition count** (per replicate, after excluding failed/rain partitions in your instrument software)

From these two numbers, the app calculates the average copies per partition (λ)
using the Poisson distribution:

```
p = positive / total
λ = -ln(1 - p)
```

A 95% confidence interval is calculated using the standard delta-method normal
approximation for Poisson-derived λ.

**Normalization & Copy Number:**
```
NF (normalization factor) = geometric mean of reference-locus λ values
Ratio = λ(target) / NF
CN (copy number) = ploidy(reference) × Ratio
Fold Change = Ratio(sample) / Ratio(control)
```

Unlike qPCR, dPCR is an **absolute, endpoint** measurement — no amplification
efficiency correction (Pfaffl) or standard curve is required, since efficiency
does not affect the binary positive/negative call at reaction endpoint.

**Quality control:** Replicates with total accepted partitions below your
configured threshold are flagged (dMIQE recommends ≥10,000 for ddPCR
instruments as a general guideline — check your instrument's own
specifications).

**Statistical testing:** Because λ and Ratio are already on a linear scale (no
log transform needed, unlike ΔCt), Shapiro-Wilk → Levene → t-test/Welch/Mann-Whitney
(2 groups) or ANOVA/Kruskal-Wallis + post-hoc (≥3 groups) are applied directly
to the Ratio values.

**Reference:** Huggett JF et al. *The digital MIQE guidelines: Minimum Information
for Publication of Quantitative Digital PCR Experiments.* Clin Chem 2013;
updated Huggett et al. 2020.

**Disclaimer:** For research and educational use only. Not validated for
clinical diagnostic decision-making.

**Contact:** mailtoburhanettin@gmail.com
""")
        else:
            st.markdown("""
### dPCR analizi qPCR'den (GeneQuantify) nasıl farklıdır?

Dijital PCR, her örneği binlerce mikro-reaksiyona (damlacık veya kuyucuk) böler.
Her partisyon **pozitif** (hedef var) veya **negatif** (hedef yok) olarak sonuçlanır.
Cq değeri yerine buraya girilen veri:

- **Pozitif partisyon sayısı** (replikat başına)
- **Toplam kabul edilen partisyon sayısı** (replikat başına; cihaz yazılımınızda
  başarısız/rain partisyonlar hariç tutulduktan sonra)

Bu iki sayıdan, Poisson dağılımı kullanılarak partisyon başına ortalama kopya
sayısı (λ) hesaplanır:

```
p = pozitif / toplam
λ = -ln(1 - p)
```

%95 güven aralığı, Poisson-türevli λ için standart delta-method normal
yaklaşımı kullanılarak hesaplanır.

**Normalizasyon & Kopya Sayısı:**
```
NF (normalizasyon faktörü) = referans lokus λ değerlerinin geometrik ortalaması
Oran = λ(hedef) / NF
CN (kopya sayısı) = ploidi(referans) × Oran
Kat Değişimi = Oran(örnek) / Oran(kontrol)
```

qPCR'nin aksine dPCR, **mutlak ve endpoint** bir ölçümdür — amplifikasyon
verimliliği düzeltmesi (Pfaffl) veya standart eğri gerekmez, çünkü verimlilik
reaksiyon sonundaki ikili pozitif/negatif sonucu etkilemez.

**Kalite kontrol:** Toplam kabul edilen partisyonu ayarladığınız eşiğin altında
kalan replikatlar işaretlenir (dMIQE genel öneri olarak ddPCR cihazları için
tipik olarak ≥10.000 önerir — kendi cihazınızın spesifikasyonlarını kontrol edin).

**İstatistiksel test:** λ ve Oran zaten doğrusal ölçekte olduğundan (ΔCt'nin
aksine log dönüşümüne gerek yoktur), Shapiro-Wilk → Levene → t-testi/Welch/
Mann-Whitney (2 grup) veya ANOVA/Kruskal-Wallis + post-hoc (≥3 grup) doğrudan
Oran değerlerine uygulanır.

**Referans:** Huggett JF ve ark. *The digital MIQE guidelines.* Clin Chem 2013;
güncelleme Huggett ve ark. 2020.

**Sorumluluk Reddi:** Yalnızca araştırma ve eğitim amaçlıdır. Klinik tanı
kararları için doğrulanmamıştır.

**İletişim:** mailtoburhanettin@gmail.com
""")
    show_guide()

st.sidebar.divider()
st.sidebar.link_button(
    _t['sidebar_github_btn'],
    "https://github.com/burhanettiny/GeneQuantify",
    use_container_width=True
)
st.sidebar.link_button(
    _t['sidebar_sister_tool'],
    "https://GeneQuantify.streamlit.app/",
    use_container_width=True
)
st.sidebar.caption("AbsoluteGene — GPL-3.0 | mailtoburhanettin@gmail.com")

tab_data, tab_results, tab_report = st.tabs([
    f"📥 {_t['tab_data']}", f"📊 {_t['tab_results']}", f"📄 {_t['tab_report']}"
])

# ═══════════════════════════════════════════════════════════════════════════════
# DATA CONTAINERS
# ═══════════════════════════════════════════════════════════════════════════════
input_values_table = []   # per-replicate rows (for display + PDF)
data = []                 # per gene/group summary rows (ratio, CN, FC)
stats_data = []            # pairwise stats rows

with tab_data:
    # ── Study Design ──────────────────────────────────────────────────────────
    with st.container(border=True):
        st.markdown(f"**{_t['study_design']}**")
        sd_c1, sd_c2, sd_c3 = st.columns(3)
        with sd_c1:
            num_target_genes = st.number_input(_t['num_target_genes'], min_value=1, step=1, key="gene_count")
        with sd_c2:
            num_patient_groups = st.number_input(_t['num_patient_groups'], min_value=1, step=1, key="patient_count")
        with sd_c3:
            num_ref_genes = st.number_input(_t['num_ref_genes'], min_value=1, max_value=10, value=1, step=1,
                                             key="num_ref_genes", help=_t['ref_gene_help'])
        sd_c4, sd_c5, sd_c6 = st.columns(3)
        with sd_c4:
            ploidy = st.number_input(_t['ploidy_label'], min_value=1, max_value=10, value=2, step=1,
                                      key="ploidy", help=_t['ploidy_help'])
        with sd_c5:
            partition_vol_nl = st.number_input(_t['partition_vol_label'], min_value=0.01, max_value=100.0,
                                                value=0.85, step=0.01, format="%.2f",
                                                key="partition_vol", help=_t['partition_vol_help'])
        with sd_c6:
            qc_min_partitions = st.number_input(_t['qc_min_partitions'], min_value=100, max_value=100000,
                                                 value=10000, step=500, key="qc_min",
                                                 help=_t['qc_min_partitions_help'])

    # ── Outlier Detection Settings ────────────────────────────────────────────
    with st.container(border=True):
        st.markdown(_t['outlier_section_title'])
        out_c1, out_c2 = st.columns([1, 2])
        with out_c1:
            outlier_enabled = st.checkbox(_t['outlier_enable'], value=True, key="outlier_enabled",
                                           help=_t['outlier_enable_help'])
            outlier_method = st.radio(_t['outlier_method_label'], options=["Grubbs", "IQR"],
                                       key="outlier_method", help=_t['outlier_method_help'])
        with out_c2:
            if outlier_method == "Grubbs":
                grubbs_alpha = st.number_input(_t['outlier_alpha_label'], min_value=0.01, max_value=0.10,
                                                value=0.05, step=0.01, format="%.2f", key="grubbs_alpha")
                iqr_multiplier = 1.5
            else:
                iqr_multiplier = st.number_input(_t['outlier_iqr_label'], min_value=1.0, max_value=3.0,
                                                  value=1.5, step=0.25, format="%.2f", key="iqr_mult")
                grubbs_alpha = 0.05

    st.divider()
    st.markdown(
        f"<div style='font-size:15px;font-weight:700;color:#004d40;margin-bottom:6px;'>"
        f"{_t['patient_data_header']}</div>", unsafe_allow_html=True
    )
    st.caption(_t['input_format_info'])

    # ── Helper: sync target + all reference loci, compute lambda/NF/ratio ──────
    def sync_and_compute(target_pos_txt, target_tot_txt, ref_pos_txts, ref_tot_txts, label, key_prefix):
        """
        Parses target & reference-locus positive/total partition text areas,
        aligns them to a common replicate count, computes Poisson lambda for
        each locus/replicate, applies QC (min partitions) and outlier
        detection (cascading across target + all reference loci), then
        computes the normalization factor (geometric mean of reference
        lambdas) and the normalized Ratio = lambda(target) / NF per kept
        replicate.
        Returns None if any required field is empty.
        """
        pos_t = parse_input_data(target_pos_txt)
        tot_t = parse_input_data(target_tot_txt)
        ref_pos = [parse_input_data(x) for x in ref_pos_txts]
        ref_tot = [parse_input_data(x) for x in ref_tot_txts]

        lengths = [len(pos_t), len(tot_t)] + [len(a) for a in ref_pos] + [len(a) for a in ref_tot]
        if any(l == 0 for l in lengths):
            return None

        n_common = min(lengths)
        if len(set(lengths)) > 1:
            st.warning(f"⚠️ {label}: " + (
                f"Farklı replikat sayıları tespit edildi — analiz n={n_common} ile devam edecek."
                if language_code == "tr" else
                f"Unequal replicate counts detected — analysis will proceed with n={n_common}."
            ))

        pos_t, tot_t = pos_t[:n_common], tot_t[:n_common]
        ref_pos = [a[:n_common] for a in ref_pos]
        ref_tot = [a[:n_common] for a in ref_tot]

        # --- target lambda/status ---
        lam_t, status_t = [], []
        for p, t in zip(pos_t, tot_t):
            lam, cl, ch, s = poisson_lambda(p, t)
            lam_t.append(lam if lam is not None else np.nan)
            status_t.append(s)
        lam_t = np.array(lam_t)
        qc_t = tot_t < qc_min_partitions

        # --- reference loci lambda/status ---
        lam_refs, status_refs, qc_refs = [], [], []
        for r in range(len(ref_pos)):
            lam_r, status_r = [], []
            for p, t in zip(ref_pos[r], ref_tot[r]):
                lam, cl, ch, s = poisson_lambda(p, t)
                lam_r.append(lam if lam is not None else np.nan)
                status_r.append(s)
            lam_refs.append(np.array(lam_r))
            status_refs.append(status_r)
            qc_refs.append(ref_tot[r] < qc_min_partitions)

        # --- combined validity mask (all loci must be ok + QC-pass at a replicate) ---
        valid_mask = np.array([s == "ok" for s in status_t]) & (~qc_t)
        for r in range(len(ref_pos)):
            valid_mask &= np.array([s == "ok" for s in status_refs[r]]) & (~qc_refs[r])

        n_qc_fail = int(np.sum(qc_t) + sum(int(np.sum(q)) for q in qc_refs))
        n_saturated = int(sum(1 for s in status_t if s == "saturated") +
                          sum(sum(1 for s in sr if s == "saturated") for sr in status_refs))
        if n_qc_fail > 0:
            st.warning(_t['qc_fail_warning'].format(n=n_qc_fail, thr=int(qc_min_partitions)))
        if n_saturated > 0:
            st.warning(_t['saturation_warning'].format(n=n_saturated))

        # --- cascading outlier detection: target first, then each ref locus ---
        excluded_idx = set()
        valid_positions = [i for i in range(n_common) if valid_mask[i]]

        if outlier_enabled and len(valid_positions) >= 3:
            target_vals = lam_t[valid_positions]
            _, excl_local = render_outlier_ui(
                target_vals, f"{label} — {_t['target_gene']}", f"{key_prefix}_tgt",
                outlier_method, grubbs_alpha, iqr_multiplier
            )
            excluded_idx.update(valid_positions[i] for i in excl_local)

        valid_positions = [i for i in valid_positions if i not in excluded_idx]
        for r in range(len(ref_pos)):
            if outlier_enabled and len(valid_positions) >= 3:
                ref_vals = lam_refs[r][valid_positions]
                _, excl_local = render_outlier_ui(
                    ref_vals, f"{label} — {_t['reference_gene']} {r+1}", f"{key_prefix}_ref{r}",
                    outlier_method, grubbs_alpha, iqr_multiplier
                )
                new_excl = [valid_positions[i] for i in excl_local]
                excluded_idx.update(new_excl)
                valid_positions = [i for i in valid_positions if i not in new_excl]

        kept = sorted(valid_positions)

        if kept:
            ref_matrix_kept = np.vstack([lam_refs[r][kept] for r in range(len(ref_pos))])
            NF_kept = np.array([geometric_mean(ref_matrix_kept[:, k]) for k in range(len(kept))])
            lam_t_kept = lam_t[kept]
            ratio_kept = lam_t_kept / NF_kept
        else:
            ref_matrix_kept = np.zeros((len(ref_pos), 0))
            NF_kept = np.array([])
            lam_t_kept = np.array([])
            ratio_kept = np.array([])

        detail_rows = []
        for i in range(n_common):
            detail_rows.append({
                "replicate": i + 1,
                "target_pos": pos_t[i], "target_tot": tot_t[i],
                "target_lambda": lam_t[i], "target_status": status_t[i],
                "ref_lambdas": [lam_refs[r][i] for r in range(len(ref_pos))],
                "qc_fail": bool(qc_t[i]) or any(bool(qc_refs[r][i]) for r in range(len(ref_pos))),
                "excluded": i in excluded_idx,
                "used": i in kept,
            })

        return {
            "n_common": n_common, "kept": kept,
            "lam_target_kept": lam_t_kept, "NF_kept": NF_kept, "ratio_kept": ratio_kept,
            "ref_matrix_kept": ref_matrix_kept,
            "detail_rows": detail_rows,
        }

    # ── Main per-gene loop ────────────────────────────────────────────────────
    for i in range(num_target_genes):
        st.markdown(
            f"<h4 style='margin-top:14px;margin-bottom:4px;color:#004d40;'>"
            f"🧬 {_t['target_gene']} {i+1}</h4>", unsafe_allow_html=True
        )

        # ── Control group input ──────────────────────────────────────────────
        st.markdown(f"**{_t['control_group']} — {_t['target_gene']} {i+1}**")
        cc1, cc2 = st.columns(2)
        with cc1:
            ctrl_target_pos_txt = st.text_area(
                f"Control {i+1} — {_t['positive_partitions']} ({_t['target_gene']})",
                value=st.session_state.get(f"ctrl_tgt_pos_{i}", ""), key=f"ctrl_tgt_pos_{i}"
            )
        with cc2:
            ctrl_target_tot_txt = st.text_area(
                f"Control {i+1} — {_t['total_partitions']} ({_t['target_gene']})",
                value=st.session_state.get(f"ctrl_tgt_tot_{i}", ""), key=f"ctrl_tgt_tot_{i}"
            )

        ctrl_ref_pos_txts, ctrl_ref_tot_txts = [], []
        for r in range(num_ref_genes):
            rc1, rc2 = st.columns(2)
            ref_lbl = f"{_t['reference_gene']} {r+1}" if num_ref_genes > 1 else _t['reference_gene']
            with rc1:
                rp = st.text_area(f"Control {i+1} — {_t['positive_partitions']} ({ref_lbl})",
                                   value=st.session_state.get(f"ctrl_ref_pos_{i}_{r}", ""),
                                   key=f"ctrl_ref_pos_{i}_{r}")
            with rc2:
                rt = st.text_area(f"Control {i+1} — {_t['total_partitions']} ({ref_lbl})",
                                   value=st.session_state.get(f"ctrl_ref_tot_{i}_{r}", ""),
                                   key=f"ctrl_ref_tot_{i}_{r}")
            ctrl_ref_pos_txts.append(rp)
            ctrl_ref_tot_txts.append(rt)

        ctrl_result = sync_and_compute(
            ctrl_target_pos_txt, ctrl_target_tot_txt, ctrl_ref_pos_txts, ctrl_ref_tot_txts,
            f"{_t['control_group']} {i+1}", f"ctrl_{i}"
        )
        if ctrl_result is None or len(ctrl_result["kept"]) == 0:
            st.error(_t['warning_empty_input'])
            continue

        # ── geNorm-style reference stability (control) ──────────────────────
        if num_ref_genes >= 2 and ctrl_result["ref_matrix_kept"].shape[1] >= 2:
            m_vals = compute_stability_m(ctrl_result["ref_matrix_kept"])
            st.markdown(f"##### 📊 {_t['genorm_title']} — {_t['control_group']} {i+1}")
            stab_cols = st.columns(num_ref_genes)
            for r, col in enumerate(stab_cols):
                with col:
                    st.metric(label=f"{_t['reference_gene']} {r+1}", value=f"M = {m_vals[r]:.3f}")
                    if m_vals[r] < 0.5:
                        st.caption(f"✅ {_t['stable']}")
                    elif m_vals[r] < 1.0:
                        st.caption(f"⚠️ {_t['borderline']}")
                    else:
                        st.caption(f"❌ {_t['unstable']}")

        # ── Store control replicate rows in input_values_table ───────────────
        for row in ctrl_result["detail_rows"]:
            status_label = (
                f"{_t['outlier_excluded_yes']} ({outlier_method})" if row["excluded"]
                else (_t['qc_fail'] if row["qc_fail"] else
                      (_t['qc_saturated'] if row["target_status"] == "saturated" else _t['outlier_excluded_no']))
            )
            input_values_table.append({
                "__gene__": f"Gene {i+1}", "Grup": "Control",
                "__replicate__": row["replicate"],
                "__positive__": row["target_pos"], "__total__": row["target_tot"],
                "__lambda__": round(row["target_lambda"], 5) if not np.isnan(row["target_lambda"]) else None,
                "Outlier Excluded": status_label,
                "__used__": row["used"],
            })

        avg_ctrl_ratio = float(np.mean(ctrl_result["ratio_kept"])) if len(ctrl_result["ratio_kept"]) > 0 else None

        # ── Patient groups ────────────────────────────────────────────────────
        for j in range(num_patient_groups):
            st.markdown(f"**{_t['patient_group']} {j+1} — {_t['target_gene']} {i+1}**")
            pc1, pc2 = st.columns(2)
            with pc1:
                smp_target_pos_txt = st.text_area(
                    f"Group {j+1} — {_t['positive_partitions']} ({_t['target_gene']} {i+1})",
                    value=st.session_state.get(f"smp_tgt_pos_{i}_{j}", ""), key=f"smp_tgt_pos_{i}_{j}"
                )
            with pc2:
                smp_target_tot_txt = st.text_area(
                    f"Group {j+1} — {_t['total_partitions']} ({_t['target_gene']} {i+1})",
                    value=st.session_state.get(f"smp_tgt_tot_{i}_{j}", ""), key=f"smp_tgt_tot_{i}_{j}"
                )

            smp_ref_pos_txts, smp_ref_tot_txts = [], []
            for r in range(num_ref_genes):
                rc1, rc2 = st.columns(2)
                ref_lbl = f"{_t['reference_gene']} {r+1}" if num_ref_genes > 1 else _t['reference_gene']
                with rc1:
                    rp = st.text_area(f"Group {j+1} — {_t['positive_partitions']} ({ref_lbl})",
                                       value=st.session_state.get(f"smp_ref_pos_{i}_{j}_{r}", ""),
                                       key=f"smp_ref_pos_{i}_{j}_{r}")
                with rc2:
                    rt = st.text_area(f"Group {j+1} — {_t['total_partitions']} ({ref_lbl})",
                                       value=st.session_state.get(f"smp_ref_tot_{i}_{j}_{r}", ""),
                                       key=f"smp_ref_tot_{i}_{j}_{r}")
                smp_ref_pos_txts.append(rp)
                smp_ref_tot_txts.append(rt)

            smp_result = sync_and_compute(
                smp_target_pos_txt, smp_target_tot_txt, smp_ref_pos_txts, smp_ref_tot_txts,
                f"{_t['patient_group']} {j+1} — Gene {i+1}", f"smp_{i}_{j}"
            )
            if smp_result is None or len(smp_result["kept"]) == 0:
                st.error(_t['warning_empty_input'])
                continue

            if num_ref_genes >= 2 and smp_result["ref_matrix_kept"].shape[1] >= 2:
                m_vals_s = compute_stability_m(smp_result["ref_matrix_kept"])
                st.markdown(f"##### 📊 {_t['genorm_title']} — {_t['patient_group']} {j+1}")
                stab_cols_s = st.columns(num_ref_genes)
                for r, col in enumerate(stab_cols_s):
                    with col:
                        st.metric(label=f"{_t['reference_gene']} {r+1}", value=f"M = {m_vals_s[r]:.3f}")
                        if m_vals_s[r] < 0.5:
                            st.caption(f"✅ {_t['stable']}")
                        elif m_vals_s[r] < 1.0:
                            st.caption(f"⚠️ {_t['borderline']}")
                        else:
                            st.caption(f"❌ {_t['unstable']}")

            for row in smp_result["detail_rows"]:
                status_label = (
                    f"{_t['outlier_excluded_yes']} ({outlier_method})" if row["excluded"]
                    else (_t['qc_fail'] if row["qc_fail"] else
                          (_t['qc_saturated'] if row["target_status"] == "saturated" else _t['outlier_excluded_no']))
                )
                input_values_table.append({
                    "__gene__": f"Gene {i+1}", "Grup": f"Group {j+1}",
                    "__replicate__": row["replicate"],
                    "__positive__": row["target_pos"], "__total__": row["target_tot"],
                    "__lambda__": round(row["target_lambda"], 5) if not np.isnan(row["target_lambda"]) else None,
                    "Outlier Excluded": status_label,
                    "__used__": row["used"],
                })

            avg_smp_ratio = float(np.mean(smp_result["ratio_kept"]))
            fold_change = avg_smp_ratio / avg_ctrl_ratio if avg_ctrl_ratio else float("nan")
            cn_ctrl = ploidy * avg_ctrl_ratio if avg_ctrl_ratio else float("nan")
            cn_smp = ploidy * avg_smp_ratio

            if fold_change >= 1.5:
                regulation = _t['upregulated']
            elif fold_change <= 0.67:
                regulation = _t['downregulated']
            else:
                regulation = _t['no_change']

            st.markdown(f"#### {_t['method_comparison']} — {_t['target_gene']} {i+1} / {_t['patient_group']} {j+1}")
            rcol1, rcol2, rcol3 = st.columns(3)
            rcol1.metric(_t['ratio_col'], f"{avg_smp_ratio:.4f}")
            rcol2.metric(_t['cn_col'], f"{cn_smp:.3f}")
            rcol3.metric(_t['fc_col'], f"{fold_change:.4f}", delta=regulation)

            # ── Statistics (Control vs this patient group), directly on Ratio ──
            control_ratios = ctrl_result["ratio_kept"]
            sample_ratios = smp_result["ratio_kept"]
            n_ctrl, n_smp = len(control_ratios), len(sample_ratios)

            if n_ctrl < 2 or n_smp < 2:
                test_pvalue, test_method, test_type, significance = float('nan'), "N/A (n < 2)", "—", "—"
                control_normal = sample_normal = True
                equal_variance = True
                shapiro_control = shapiro_sample = type('SW', (), {'statistic': float('nan'), 'pvalue': float('nan')})()
                levene_test = type('LV', (), {'statistic': float('nan'), 'pvalue': float('nan')})()
            else:
                _MIN_N_SHAPIRO = 8
                if n_ctrl >= _MIN_N_SHAPIRO and n_smp >= _MIN_N_SHAPIRO:
                    shapiro_control = stats.shapiro(control_ratios)
                    shapiro_sample = stats.shapiro(sample_ratios)
                    control_normal = shapiro_control.pvalue > 0.05
                    sample_normal = shapiro_sample.pvalue > 0.05
                else:
                    shapiro_control = shapiro_sample = type('SW', (), {'statistic': float('nan'), 'pvalue': float('nan')})()
                    control_normal = sample_normal = True

                try:
                    levene_test = stats.levene(control_ratios, sample_ratios)
                    equal_variance = (levene_test.pvalue > 0.05) if not np.isnan(levene_test.pvalue) else True
                except Exception:
                    levene_test = type('LV', (), {'statistic': float('nan'), 'pvalue': float('nan')})()
                    equal_variance = True

                try:
                    if control_normal and sample_normal:
                        if equal_variance:
                            test_pvalue = stats.ttest_ind(control_ratios, sample_ratios).pvalue
                            test_method = _t['t_test']
                        else:
                            test_pvalue = stats.ttest_ind(control_ratios, sample_ratios, equal_var=False).pvalue
                            test_method = _t['welch_t_test']
                        test_type = _t['parametric']
                    else:
                        test_pvalue = stats.mannwhitneyu(control_ratios, sample_ratios, alternative='two-sided').pvalue
                        test_method = _t['mann_whitney_u_test']
                        test_type = _t['non_parametric']
                except Exception:
                    test_pvalue, test_method, test_type = float('nan'), "Error", "—"

                significance = (_t['significant'] if (not np.isnan(test_pvalue) and test_pvalue < 0.05)
                                 else (_t['insignificant'] if not np.isnan(test_pvalue) else "—"))

            with st.expander(f"{_t['stat_decision_title']} — {_t['target_gene']} {i+1} / Group {j+1}", expanded=False):
                st.markdown(_t['stat_decision_steps'])
                if n_ctrl >= 8 and n_smp >= 8:
                    sw_ctrl_sym = "✅" if control_normal else "❌"
                    sw_smp_sym = "✅" if sample_normal else "❌"
                    st.markdown(
                        f"{_t['stat_shapiro_title']}  \n"
                        f"- Control: W={shapiro_control.statistic:.4f}, p={shapiro_control.pvalue:.4f} {sw_ctrl_sym} "
                        f"{_t['stat_normal'] if control_normal else _t['stat_nonnormal']}  \n"
                        f"- Group {j+1}: W={shapiro_sample.statistic:.4f}, p={shapiro_sample.pvalue:.4f} {sw_smp_sym} "
                        f"{_t['stat_normal'] if sample_normal else _t['stat_nonnormal']}"
                    )
                else:
                    st.info(f"ℹ️ Shapiro-Wilk skipped (n={min(n_ctrl, n_smp)} < 8). Normality assumed.")

                if control_normal and sample_normal:
                    lev_sym = "✅" if equal_variance else "⚠️"
                    st.markdown(
                        f"{_t['stat_levene_title']}  \n"
                        f"- F={levene_test.statistic:.4f}, p={levene_test.pvalue:.4f} {lev_sym} "
                        f"{_t['stat_equal_var'] if equal_variance else _t['stat_unequal_var']}"
                    )
                else:
                    st.markdown(_t['stat_levene_skipped'])

                if not control_normal or not sample_normal:
                    reason = _t['stat_reason_nonnormal']
                elif equal_variance:
                    reason = _t['stat_reason_normal_equal']
                else:
                    reason = _t['stat_reason_normal_unequal']

                st.success(
                    f"{_t['stat_selected_test']} {test_method}  \n"
                    f"{_t['stat_reason']} {reason}  \n"
                    f"{_t['stat_result']} p = {test_pvalue:.4f} → **{significance}**"
                )
                if num_patient_groups >= 2:
                    st.caption(_t['stat_multigroup_note'])

            stats_data.append({
                "__gene__": f"Gene {i+1}", "__group__": f"Group {j+1}",
                "__test_type__": test_type, "__test_method__": test_method,
                "__pvalue__": test_pvalue, "__significance__": significance,
                "Comparison": f"Control vs Group {j+1}",
                "__ratio_ctrl__": control_ratios, "__ratio_smp__": sample_ratios,
            })

            data.append({
                "__gene__": f"Gene {i+1}", "__group__": f"Group {j+1}",
                "__ratio_ctrl__": avg_ctrl_ratio, "__ratio_smp__": avg_smp_ratio,
                "__cn_ctrl__": cn_ctrl, "__cn_smp__": cn_smp,
                "__fc__": fold_change, "__regulation__": regulation,
                "__n_ctrl__": n_ctrl, "__n_smp__": n_smp,
            })

# ═══════════════════════════════════════════════════════════════════════════════
# MULTI-GROUP ANALYSIS (≥3 groups per gene) — omnibus + post-hoc + correction
# ═══════════════════════════════════════════════════════════════════════════════
multigroup_results = []

for i in range(num_target_genes):
    gene_label = f"Gene {i+1}"
    gene_rows = [r for r in stats_data if r["__gene__"] == gene_label]
    if not gene_rows:
        continue

    ctrl_ratios = gene_rows[0]["__ratio_ctrl__"]
    group_ratio_map = {r["__group__"]: r["__ratio_smp__"] for r in gene_rows}

    all_groups = [list(ctrl_ratios)] + [list(v) for v in group_ratio_map.values()]
    all_group_names = ["Control"] + list(group_ratio_map.keys())
    n_groups = len(all_groups)

    if n_groups < 3:
        multigroup_results.append({
            "gene": gene_label, "n_groups": n_groups, "note": "2-group",
            "omnibus_test": "—", "omnibus_p": None, "posthoc_rows": []
        })
        continue

    normality_ok = all((len(g) < 8 or stats.shapiro(g).pvalue > 0.05) for g in all_groups if len(g) >= 3)
    levene_p = stats.levene(*all_groups).pvalue if n_groups >= 2 else 1.0
    variance_ok = levene_p > 0.05

    if normality_ok and variance_ok:
        omnibus_stat, omnibus_p = stats.f_oneway(*all_groups)
        omnibus_test, omnibus_type, posthoc_method = "One-way ANOVA", "parametric", "Tukey HSD"
    elif normality_ok and not variance_ok:
        try:
            from scipy.stats import alexandergovern
            result = alexandergovern(*all_groups)
            omnibus_p, omnibus_stat = result.pvalue, result.statistic
        except Exception:
            omnibus_stat, omnibus_p = stats.f_oneway(*all_groups)
        omnibus_test, omnibus_type, posthoc_method = "Welch ANOVA", "parametric", "Games-Howell (approx.)"
    else:
        omnibus_stat, omnibus_p = stats.kruskal(*all_groups)
        omnibus_test, omnibus_type, posthoc_method = "Kruskal-Wallis", "non-parametric", "Dunn (Mann-Whitney U)"

    omnibus_sig = _t['multigroup_significant'] if omnibus_p < 0.05 else _t['multigroup_not_significant']

    pairs, raw_pvals = [], []
    for a in range(n_groups):
        for b in range(a + 1, n_groups):
            g_a, g_b = all_groups[a], all_groups[b]
            if omnibus_type == "parametric" and variance_ok:
                p = stats.ttest_ind(g_a, g_b).pvalue
            elif omnibus_type == "parametric" and not variance_ok:
                p = stats.ttest_ind(g_a, g_b, equal_var=False).pvalue
            else:
                p = stats.mannwhitneyu(g_a, g_b, alternative="two-sided").pvalue
            pairs.append((all_group_names[a], all_group_names[b]))
            raw_pvals.append(p)

    n_tests = len(raw_pvals)
    bonf_pvals = [min(p * n_tests, 1.0) for p in raw_pvals]
    ranked = sorted(range(n_tests), key=lambda k: raw_pvals[k])
    fdr_pvals = [1.0] * n_tests
    for rank, idx in enumerate(ranked):
        fdr_pvals[idx] = min(raw_pvals[idx] * n_tests / (rank + 1), 1.0)
    for k in range(n_tests - 2, -1, -1):
        fdr_pvals[ranked[k]] = min(fdr_pvals[ranked[k]], fdr_pvals[ranked[k + 1]])

    posthoc_rows = []
    for idx, (pa, pb) in enumerate(pairs):
        posthoc_rows.append({
            "Comparison": f"{pa} vs {pb}", "Raw p": round(raw_pvals[idx], 4),
            "Bonferroni p": round(bonf_pvals[idx], 4), "FDR p (B-H)": round(fdr_pvals[idx], 4),
            "Sig (raw)": "✅" if raw_pvals[idx] < 0.05 else "—",
            "Sig (Bonferroni)": "✅" if bonf_pvals[idx] < 0.05 else "—",
            "Sig (FDR)": "✅" if fdr_pvals[idx] < 0.05 else "—",
        })

    multigroup_results.append({
        "gene": gene_label, "n_groups": n_groups,
        "omnibus_test": omnibus_test, "omnibus_type": omnibus_type, "omnibus_p": omnibus_p,
        "omnibus_sig": omnibus_sig, "posthoc_method": posthoc_method, "posthoc_rows": posthoc_rows,
        "normality_ok": normality_ok, "variance_ok": variance_ok, "note": None
    })

# ═══════════════════════════════════════════════════════════════════════════════
# RESULTS TAB
# ═══════════════════════════════════════════════════════════════════════════════
with tab_results:

    # ── Multi-group display ───────────────────────────────────────────────────
    if any(r["n_groups"] >= 3 for r in multigroup_results):
        st.markdown("---")
        st.markdown(_t['multigroup_title'])
        for res in multigroup_results:
            if res["n_groups"] < 3:
                continue
            st.markdown(f"### 🧬 {res['gene']} — {res['n_groups']} {_t['patient_group'].replace('🩸 ', '')}")
            if res["normality_ok"] and res["variance_ok"]:
                st.success(_t['multigroup_decision_normal_equal'])
            elif res["normality_ok"] and not res["variance_ok"]:
                st.warning(_t['multigroup_decision_normal_unequal'])
            else:
                st.warning(_t['multigroup_decision_nonnormal'])

            omni_col1, omni_col2, omni_col3 = st.columns(3)
            omni_col1.metric(_t['multigroup_omnibus_test'], res["omnibus_test"])
            omni_col2.metric(_t['multigroup_pvalue'], f"{res['omnibus_p']:.4f}")
            omni_col3.metric(_t['multigroup_result'], res["omnibus_sig"])
            if res["omnibus_p"] >= 0.05:
                st.info(_t['multigroup_omnibus_ns'])

            st.markdown(f"{_t['multigroup_posthoc_label']} {res['posthoc_method']}")
            ph_df = pd.DataFrame(res["posthoc_rows"])
            st.dataframe(ph_df, use_container_width=True)

            fig_ph = go.Figure()
            comparisons = [r["Comparison"] for r in res["posthoc_rows"]]
            fig_ph.add_trace(go.Bar(name="Raw p", x=comparisons, y=[r["Raw p"] for r in res["posthoc_rows"]], marker_color="#00796b"))
            fig_ph.add_trace(go.Bar(name="Bonferroni p", x=comparisons, y=[r["Bonferroni p"] for r in res["posthoc_rows"]], marker_color="#f57c00"))
            fig_ph.add_trace(go.Bar(name="FDR p (B-H)", x=comparisons, y=[r["FDR p (B-H)"] for r in res["posthoc_rows"]], marker_color="#455a64"))
            fig_ph.add_hline(y=0.05, line_dash="dash", line_color="red", annotation_text="α = 0.05")
            fig_ph.update_layout(barmode="group", title=f"{res['gene']} — Post-hoc p-values", yaxis_title="p-value", height=350)
            st.plotly_chart(fig_ph, use_container_width=True, key=f"posthoc_{res['gene']}")

            ph_csv = ph_df.to_csv(index=False).encode("utf-8")
            st.download_button(f"{_t['multigroup_dl_button']} {res['gene']}", data=ph_csv,
                                file_name=f"posthoc_{res['gene'].replace(' ', '_')}.csv", mime="text/csv",
                                key=f"ph_dl_{res['gene']}")
    elif num_patient_groups >= 2 and multigroup_results:
        st.markdown("---")
        st.info(_t['multigroup_2group_note'])

    # ── Input data table ──────────────────────────────────────────────────────
    if input_values_table:
        st.subheader(_t['gr_tbl'])
        _rename = {
            "__gene__": _t['target_gene'], "Grup": "Group", "__replicate__": _t['sample_number'],
            "__positive__": _t['positive_partitions'], "__total__": _t['total_partitions'],
            "__lambda__": _t['lambda_col'], "Outlier Excluded": "Status",
        }
        display_df = pd.DataFrame(input_values_table).drop(columns=["__used__"]).rename(columns=_rename)
        st.dataframe(display_df, use_container_width=True)
        csv = display_df.to_csv(index=False).encode("utf-8")
        st.download_button(_t['download_csv'], data=csv, file_name="dpcr_input_data.csv", mime="text/csv", key="dl_input_csv")

    # ── Results summary table ─────────────────────────────────────────────────
    if data:
        st.subheader(_t['nil_mine'])
        _rename2 = {
            "__gene__": _t['target_gene'], "__group__": "Group",
            "__ratio_ctrl__": f"{_t['ratio_col']} (Control)", "__ratio_smp__": f"{_t['ratio_col']} (Sample)",
            "__cn_ctrl__": f"{_t['cn_col']} (Control)", "__cn_smp__": f"{_t['cn_col']} (Sample)",
            "__fc__": _t['fc_col'], "__regulation__": _t['regulation_status'],
            "__n_ctrl__": "n (Control)", "__n_smp__": "n (Sample)",
        }
        res_df = pd.DataFrame(data).rename(columns=_rename2)
        st.dataframe(res_df, use_container_width=True)
        csv2 = res_df.to_csv(index=False).encode("utf-8")
        st.download_button(_t['download_csv'], data=csv2, file_name="dpcr_results.csv", mime="text/csv", key="dl_res_csv")

    # ── Statistics table ───────────────────────────────────────────────────────
    if stats_data:
        st.subheader(_t['statistical_results'])
        _rename3 = {
            "__gene__": _t['target_gene'], "__group__": "Group", "__test_type__": _t['test_type'],
            "__test_method__": _t['test_method'], "__pvalue__": _t['test_pvalue'], "__significance__": _t['significance'],
        }
        stats_df = pd.DataFrame(stats_data).drop(columns=["__ratio_ctrl__", "__ratio_smp__"]).rename(columns=_rename3)
        st.dataframe(stats_df, use_container_width=True)
        csv3 = stats_df.to_csv(index=False).encode("utf-8")
        st.download_button(_t['download_csv'], data=csv3, file_name="dpcr_statistics.csv", mime="text/csv", key="dl_stats_csv")

    # ── Multi-gene p-value correction ─────────────────────────────────────────
    if stats_data and num_target_genes >= 2:
        st.markdown("---")
        st.markdown(_t['multigene_title'])
        correction_rows = [
            {"Gene": r["__gene__"], "Group": r["__group__"], "Raw p": r["__pvalue__"], "Test": r["__test_method__"]}
            for r in stats_data if r.get("__pvalue__") is not None and not np.isnan(r["__pvalue__"])
        ]
        if correction_rows:
            n_tests = len(correction_rows)
            raw_pvals = [r["Raw p"] for r in correction_rows]
            bonf = [min(p * n_tests, 1.0) for p in raw_pvals]
            ranked = sorted(range(n_tests), key=lambda k: raw_pvals[k])
            fdr = [1.0] * n_tests
            for rank, idx in enumerate(ranked):
                fdr[idx] = min(raw_pvals[idx] * n_tests / (rank + 1), 1.0)
            for k in range(n_tests - 2, -1, -1):
                fdr[ranked[k]] = min(fdr[ranked[k]], fdr[ranked[k + 1]])
            for idx, row in enumerate(correction_rows):
                row["Bonferroni p"] = round(bonf[idx], 4)
                row["FDR p (B-H)"] = round(fdr[idx], 4)
                row["Sig (raw)"] = "✅" if raw_pvals[idx] < 0.05 else "—"
                row["Sig (Bonferroni)"] = "✅" if bonf[idx] < 0.05 else "—"
                row["Sig (FDR)"] = "✅" if fdr[idx] < 0.05 else "—"

            corr_df = pd.DataFrame(correction_rows)
            st.dataframe(corr_df, use_container_width=True)

            n_raw_sig = sum(1 for p in raw_pvals if p < 0.05)
            n_fdr_sig = sum(1 for p in fdr if p < 0.05)
            sc1, sc2, sc3 = st.columns(3)
            sc1.metric(_t['multigene_sig_raw'], f"{n_raw_sig} / {n_tests}")
            sc2.metric(_t['multigene_sig_bonf'], f"{sum(1 for p in bonf if p < 0.05)} / {n_tests}")
            sc3.metric(_t['multigene_sig_fdr'], f"{n_fdr_sig} / {n_tests}")

            if n_raw_sig > n_fdr_sig:
                st.warning(_t['multigene_warning'].format(lost=n_raw_sig - n_fdr_sig))
            elif n_raw_sig == n_fdr_sig and n_raw_sig > 0:
                st.success(_t['multigene_success'].format(n=n_raw_sig))
            else:
                st.info(_t['multigene_no_sig'])

            corr_csv = corr_df.to_csv(index=False).encode("utf-8")
            st.download_button(_t['multigene_dl_button'], data=corr_csv, file_name="multi_gene_correction.csv",
                                mime="text/csv", key="multigene_corr_dl")
    elif stats_data and num_target_genes == 1:
        st.markdown("---")
        st.info(_t['multigene_1gene_note'])

    # ── Distribution plots ────────────────────────────────────────────────────
    st.markdown("---")
    plot_mode = st.radio(
        _t['dist_plot_mode_label'],
        options=[_t['dist_plot_ratio'], _t['dist_plot_lambda'], _t['dist_plot_fc']],
        index=0, horizontal=True, key="dist_plot_mode"
    )
    _mode_id = "RATIO" if plot_mode == _t['dist_plot_ratio'] else ("LAMBDA" if plot_mode == _t['dist_plot_lambda'] else "FC")

    for i in range(num_target_genes):
        gene_label = f"Gene {i+1}"
        gene_rows_st = [r for r in stats_data if r["__gene__"] == gene_label]
        if not gene_rows_st:
            continue
        st.subheader(f"{_t['target_gene']} {i+1} - {_t['distribution_graph']}")

        ctrl_ratios = gene_rows_st[0]["__ratio_ctrl__"]
        ctrl_lambdas = np.array([
            row["__lambda__"] for row in input_values_table
            if row["__gene__"] == gene_label and row["Grup"] == "Control" and row["__used__"] and row["__lambda__"] is not None
        ])

        if _mode_id == "RATIO":
            ctrl_vals = np.array(ctrl_ratios)
            y_label = _t['ratio_col']
        elif _mode_id == "LAMBDA":
            ctrl_vals = ctrl_lambdas
            y_label = _t['lambda_col']
        else:
            ctrl_vals = np.ones(len(ctrl_ratios))  # FC of control vs itself = 1
            y_label = _t['fc_col']

        fig = go.Figure()
        avg_ctrl = float(np.mean(ctrl_vals)) if len(ctrl_vals) else np.nan
        fig.add_trace(go.Scatter(x=[0.8, 1.2], y=[avg_ctrl, avg_ctrl], mode='lines',
                                  line=dict(color='black', width=4), name="Control avg"))
        fig.add_trace(go.Scatter(
            x=np.ones(len(ctrl_vals)) + np.random.uniform(-0.05, 0.05, len(ctrl_vals)),
            y=ctrl_vals, mode='markers', name="Control", marker=dict(color='#00796b')
        ))

        for j in range(num_patient_groups):
            row_match = [r for r in gene_rows_st if r["__group__"] == f"Group {j+1}"]
            if not row_match:
                continue
            smp_ratios = np.array(row_match[0]["__ratio_smp__"])
            if _mode_id == "RATIO":
                smp_vals = smp_ratios
            elif _mode_id == "LAMBDA":
                smp_vals = np.array([
                    row["__lambda__"] for row in input_values_table
                    if row["__gene__"] == gene_label and row["Grup"] == f"Group {j+1}" and row["__used__"] and row["__lambda__"] is not None
                ])
            else:
                smp_vals = smp_ratios / avg_ctrl if avg_ctrl else smp_ratios

            avg_smp = float(np.mean(smp_vals)) if len(smp_vals) else np.nan
            fig.add_trace(go.Scatter(x=[j + 1.8, j + 2.2], y=[avg_smp, avg_smp], mode='lines',
                                      line=dict(color='black', width=4), name=f"Group {j+1} avg"))
            fig.add_trace(go.Scatter(
                x=np.ones(len(smp_vals)) * (j + 2) + np.random.uniform(-0.05, 0.05, len(smp_vals)),
                y=smp_vals, mode='markers', name=f"Group {j+1}", marker=dict(color='#d84315')
            ))

        fig.update_layout(
            title=f"{_t['target_gene']} {i+1} — {y_label}",
            xaxis=dict(tickvals=[1] + [j + 2 for j in range(num_patient_groups)],
                       ticktext=["Control"] + [f"Group {j+1}" for j in range(num_patient_groups)],
                       title=_t['x_axis_title']),
            yaxis=dict(title=y_label), showlegend=True
        )
        st.plotly_chart(fig, use_container_width=True, key=f"dist_chart_{i}")

# ═══════════════════════════════════════════════════════════════════════════════
# PDF REPORT GENERATION
# ═══════════════════════════════════════════════════════════════════════════════
def create_pdf(results, stat_rows, input_df, lang, multigroup_results=None):
    T = translations[lang]
    buffer = BytesIO()
    doc = SimpleDocTemplate(buffer, pagesize=letter, leftMargin=50, rightMargin=50, topMargin=60, bottomMargin=50)
    elements = []
    styles = getSampleStyleSheet()

    title_style = ParagraphStyle('RT', parent=styles['Title'], fontName=PDF_FONT_BOLD, fontSize=20,
                                  textColor=colors.HexColor('#004d40'), spaceAfter=6, alignment=1)
    sub_style = ParagraphStyle('RS', parent=styles['Normal'], fontName=PDF_FONT, fontSize=10,
                                textColor=colors.HexColor('#555555'), spaceAfter=4, alignment=1)
    h1_style = ParagraphStyle('H1', parent=styles['Heading1'], fontName=PDF_FONT_BOLD, fontSize=13,
                               textColor=colors.HexColor('#004d40'), spaceBefore=16, spaceAfter=5)
    h2_style = ParagraphStyle('H2', parent=styles['Heading2'], fontName=PDF_FONT_BOLD, fontSize=11,
                               textColor=colors.HexColor('#00695c'), spaceBefore=10, spaceAfter=4)
    body_style = ParagraphStyle('BD', parent=styles['Normal'], fontName=PDF_FONT, fontSize=9, leading=13, spaceAfter=4)
    small_style = ParagraphStyle('SM', parent=styles['Normal'], fontName=PDF_FONT, fontSize=8, leading=11,
                                  textColor=colors.HexColor('#444444'))
    caption_style = ParagraphStyle('CA', parent=styles['Normal'], fontName=PDF_FONT, fontSize=8,
                                    textColor=colors.HexColor('#666666'), alignment=1, spaceAfter=6)
    warn_style = ParagraphStyle('WN', parent=styles['Normal'], fontName=PDF_FONT, fontSize=9, leading=13,
                                 backColor=colors.HexColor('#fff8e1'), borderPad=6, leftIndent=8, rightIndent=8, spaceAfter=6)

    def hr():
        return HRFlowable(width="100%", thickness=0.5, color=colors.HexColor('#cccccc'), spaceAfter=8, spaceBefore=4)

    def make_table(rows, col_widths=None, header=True):
        if not rows:
            return Spacer(1, 1)
        styled_rows = []
        for ri, row in enumerate(rows):
            styled_row = []
            for cell in row:
                if isinstance(cell, Flowable):
                    styled_row.append(cell)
                    continue
                cell_str = safe_str(cell)
                if ri == 0 and header:
                    p = Paragraph(cell_str, ParagraphStyle('TH', fontName=PDF_FONT_BOLD, fontSize=7,
                                                            textColor=colors.white, alignment=1))
                else:
                    p = Paragraph(cell_str, ParagraphStyle('TD', fontName=PDF_FONT, fontSize=7, alignment=1))
                styled_row.append(p)
            styled_rows.append(styled_row)
        tbl = Table(styled_rows, colWidths=col_widths, repeatRows=1 if header else 0)
        tbl_style = [
            ('FONTNAME', (0, 0), (-1, -1), PDF_FONT), ('ALIGN', (0, 0), (-1, -1), 'CENTER'),
            ('VALIGN', (0, 0), (-1, -1), 'MIDDLE'), ('GRID', (0, 0), (-1, -1), 0.3, colors.HexColor('#cccccc')),
            ('ROWBACKGROUNDS', (0, 1), (-1, -1), [colors.white, colors.HexColor('#e0f2f1')]),
            ('TOPPADDING', (0, 0), (-1, -1), 4), ('BOTTOMPADDING', (0, 0), (-1, -1), 4),
        ]
        if header:
            tbl_style.append(('BACKGROUND', (0, 0), (-1, 0), colors.HexColor('#004d40')))
        tbl.setStyle(TableStyle(tbl_style))
        return tbl

    import datetime
    now = datetime.datetime.now().strftime("%Y-%m-%d %H:%M")

    # ── COVER ──────────────────────────────────────────────────────────────────
    elements.append(Spacer(1, 40))
    elements.append(Paragraph(safe_str("AbsoluteGene"), title_style))
    elements.append(Paragraph(safe_str(T['pdf_report']), sub_style))
    elements.append(Spacer(1, 6))
    elements.append(Paragraph(safe_str(f"Generated: {now}"), sub_style))
    elements.append(Paragraph(safe_str("Related tool for qPCR: GeneQuantify — Yalçınkaya B (2026), Mol Cell Biochem, "
                                        "https://doi.org/10.1007/s11010-026-05621-y"), sub_style))
    elements.append(Spacer(1, 20))
    elements.append(hr())

    n_genes = len(set(r["__gene__"] for r in results)) if results else 0
    n_groups = len(set(r["__group__"] for r in results)) if results else 0
    n_replicates = len(input_df)
    n_excluded = sum(1 for _, row in input_df.iterrows() if str(row.get('Outlier Excluded', '')).startswith(('Yes', 'Evet')))
    summary_rows = [
        ["Parameter", "Value"],
        ["Target genes analyzed", str(n_genes)],
        ["Patient groups", str(n_groups)],
        ["Total replicates", str(n_replicates)],
        ["Excluded / flagged replicates", str(n_excluded)],
        ["Reference loci per gene", str(num_ref_genes)],
        ["Reference ploidy", str(ploidy)],
        ["Partition volume (nL)", str(partition_vol_nl)],
        ["QC minimum partitions", str(int(qc_min_partitions))],
        ["Calculation method", "Poisson statistics (dPCR, absolute quantification)"],
    ]
    elements.append(make_table(summary_rows, col_widths=[260, 200]))
    elements.append(Spacer(1, 14))
    elements.append(Paragraph(safe_str(
        "This report was generated automatically by AbsoluteGene. All calculations follow the "
        "digital MIQE (dMIQE) guidelines (Huggett et al., Clin Chem 2013; updated 2020). "
        "For research and educational use only — not validated for clinical diagnostic decision-making."
    ), small_style))
    elements.append(PageBreak())

    # ── SECTION 1: METHODS ────────────────────────────────────────────────────
    elements.append(Paragraph(safe_str("1. Methods and Analysis Settings"), h1_style))
    elements.append(hr())
    elements.append(Paragraph(safe_str("1.1 Poisson Quantification"), h2_style))
    elements.append(Paragraph(safe_str(
        "For each replicate, the fraction of positive partitions p = positive/total was converted to "
        "copies per partition (λ) using the Poisson distribution: λ = -ln(1-p). A 95% confidence "
        "interval was calculated using the delta-method normal approximation: "
        "Var(λ) ≈ p / (n·(1-p))."
    ), body_style))
    elements.append(Paragraph(safe_str("1.2 Normalization and Copy Number"), h2_style))
    elements.append(Paragraph(safe_str(
        f"Normalization factor (NF) = geometric mean of reference-locus λ values "
        f"({num_ref_genes} reference locus/loci used). Ratio = λ(target) / NF. "
        f"Copy number (CN) = ploidy(reference) × Ratio, using reference ploidy = {ploidy}. "
        f"Fold Change = Ratio(sample) / Ratio(control). Because dPCR is an absolute endpoint "
        f"measurement, no amplification efficiency correction is required."
    ), body_style))
    elements.append(Paragraph(safe_str("1.3 Quality Control"), h2_style))
    elements.append(Paragraph(safe_str(
        f"Replicates with fewer than {int(qc_min_partitions)} total accepted partitions were flagged "
        f"as low quality and excluded from λ-based calculations. Replicates with 100% positive "
        f"partitions (saturation) were also excluded, as λ is undefined at p=1; further dilution "
        f"is recommended for such samples."
    ), body_style))
    elements.append(Paragraph(safe_str("1.4 Outlier Detection"), h2_style))
    if outlier_enabled:
        method_txt = (f"Grubbs' test (alpha = {grubbs_alpha})" if outlier_method == "Grubbs"
                       else f"IQR method (multiplier k = {iqr_multiplier})")
        elements.append(Paragraph(safe_str(
            f"{method_txt} was applied to λ values (copies/partition) for the target locus and each "
            f"reference locus, cascading across all loci within a replicate. Flagged replicates were "
            f"confirmed for exclusion by the user; {n_excluded} replicate(s) were excluded/flagged in "
            f"this analysis."
        ), body_style))
    else:
        elements.append(Paragraph(safe_str("Outlier detection was disabled for this analysis."), body_style))
    elements.append(PageBreak())

    # ── SECTION 2: INPUT DATA ─────────────────────────────────────────────────
    elements.append(Paragraph(safe_str("2. Input Data"), h1_style))
    elements.append(hr())
    elements.append(Paragraph(safe_str(
        "Positive and total accepted partition counts entered by the user, with derived λ "
        "(copies/partition) values. Rows flagged as excluded were removed from ratio calculations."
    ), body_style))
    elements.append(Spacer(1, 6))
    if not input_df.empty:
        _pdf_rename = {
            "__gene__": "Gene", "Grup": "Group", "__replicate__": "Replicate #",
            "__positive__": "Positive", "__total__": "Total", "__lambda__": "Lambda (copies/partition)",
            "Outlier Excluded": "Status",
        }
        input_df_disp = input_df.drop(columns=["__used__"], errors="ignore").rename(columns=_pdf_rename)
        cols = input_df_disp.columns.tolist()
        page_w = letter[0] - 100
        cw = page_w / max(len(cols), 1)
        tbl_rows = [cols]
        for _, row in input_df_disp.iterrows():
            is_excl = str(row.get('Status', '')).startswith(('Yes', 'Evet'))
            row_cells = []
            for c in cols:
                v = row.get(c, '')
                cell_str = safe_str(str(v) if v is not None else '')
                style = ParagraphStyle('EX' if is_excl else 'TD', fontName=PDF_FONT, fontSize=7, alignment=1,
                                        textColor=colors.HexColor('#cc0000') if is_excl else colors.black)
                row_cells.append(Paragraph(cell_str, style))
            tbl_rows.append(row_cells)
        elements.append(make_table(tbl_rows, col_widths=[cw] * len(cols)))
    elements.append(PageBreak())

    # ── SECTION 3: RESULTS ────────────────────────────────────────────────────
    elements.append(Paragraph(safe_str("3. Results"), h1_style))
    elements.append(hr())
    res_cols = ["Gene", "Group", "Ratio (Ctrl)", "Ratio (Sample)", "CN (Ctrl)", "CN (Sample)", "Fold Change", "Regulation"]
    res_rows = [res_cols]
    for r in results:
        res_rows.append([
            r["__gene__"], r["__group__"],
            f"{r['__ratio_ctrl__']:.4f}" if r['__ratio_ctrl__'] is not None else "—",
            f"{r['__ratio_smp__']:.4f}",
            f"{r['__cn_ctrl__']:.3f}" if r['__cn_ctrl__'] is not None else "—",
            f"{r['__cn_smp__']:.3f}",
            f"{r['__fc__']:.4f}", r["__regulation__"],
        ])
    cw8 = (letter[0] - 100) / 8
    elements.append(make_table(res_rows, col_widths=[cw8] * 8))
    elements.append(Spacer(1, 8))

    if results:
        try:
            fig_fc, ax_fc = plt.subplots(figsize=(7, 3.5))
            labels_fc = [f"{r['__gene__']} /\n{r['__group__']}" for r in results]
            vals_fc = [r["__fc__"] for r in results]
            bars = ax_fc.bar(range(len(labels_fc)), vals_fc, color='#00796b', alpha=0.85)
            ax_fc.axhline(y=1, color='black', linestyle='--', linewidth=0.8, alpha=0.6)
            ax_fc.set_xticks(range(len(labels_fc)))
            ax_fc.set_xticklabels(labels_fc, fontsize=7)
            ax_fc.set_ylabel('Fold Change', fontsize=9)
            ax_fc.set_title('Fold Change by Gene/Group', fontsize=10, fontweight='bold')
            ax_fc.spines['top'].set_visible(False)
            ax_fc.spines['right'].set_visible(False)
            for bar in bars:
                ax_fc.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.02,
                           f'{bar.get_height():.2f}', ha='center', va='bottom', fontsize=6)
            plt.tight_layout()
            ib = BytesIO()
            plt.savefig(ib, format='png', dpi=150, bbox_inches='tight')
            plt.close()
            ib.seek(0)
            elements.append(RLImage(ib, width=460, height=230))
            elements.append(Paragraph(safe_str("Figure 1. Fold change by gene/group. Dashed line y=1 = no change vs control."), caption_style))
        except Exception:
            pass
    elements.append(PageBreak())

    # ── SECTION 4: STATISTICS ─────────────────────────────────────────────────
    elements.append(Paragraph(safe_str("4. Statistical Analysis"), h1_style))
    elements.append(hr())
    elements.append(Paragraph(safe_str(
        "Statistical tests were performed directly on the normalized Ratio values (linear scale, "
        "no log transform required). Test selection is automatic based on normality (Shapiro-Wilk, "
        "n≥8) and variance homogeneity (Levene). Significance threshold: p < 0.05."
    ), body_style))
    elements.append(Spacer(1, 6))

    _has_multigroup = multigroup_results is not None and any(r.get("n_groups", 0) >= 3 for r in multigroup_results)
    stat_cols = ["Gene", "Comparison", "Type", "Method", "p-value", "Significance"]
    stat_table_rows = [stat_cols]
    for sr in stat_rows:
        stat_table_rows.append([
            sr["__gene__"], sr["Comparison"], sr["__test_type__"], sr["__test_method__"],
            f"{sr['__pvalue__']:.4f}" if not np.isnan(sr['__pvalue__']) else "—", sr["__significance__"],
        ])
    cw6 = (letter[0] - 100) / 6
    elements.append(make_table(stat_table_rows, col_widths=[cw6] * 6))
    elements.append(Spacer(1, 10))

    if _has_multigroup:
        elements.append(Paragraph(safe_str("4.1 Multi-Group Comparison (≥3 Groups)"), h2_style))
        for mg in (multigroup_results or []):
            if mg.get("n_groups", 0) < 3:
                continue
            elements.append(Paragraph(safe_str(f"▶ {mg['gene']} — {mg['omnibus_test']} (p = {mg['omnibus_p']:.4f})"), body_style))
            ph_rows = mg.get("posthoc_rows", [])
            if ph_rows:
                ph_header = ["Comparison", "Raw p", "Bonferroni p", "FDR p (B-H)"]
                ph_table = [ph_header]
                for ph in ph_rows:
                    ph_table.append([ph["Comparison"], f"{ph['Raw p']:.4f}", f"{ph['Bonferroni p']:.4f}", f"{ph['FDR p (B-H)']:.4f}"])
                elements.append(make_table(ph_table, col_widths=[(letter[0] - 100) / 4] * 4))
                elements.append(Spacer(1, 8))
    elements.append(PageBreak())

    # ── SECTION 5: DISTRIBUTION PLOTS ─────────────────────────────────────────
    elements.append(Paragraph(safe_str("5. Normalized Ratio Distribution Plots"), h1_style))
    elements.append(hr())
    elements.append(Paragraph(safe_str(
        "Distribution of normalized Ratio (Target λ / reference NF) values per group. "
        "Each point = one biological replicate; horizontal bars = group means."
    ), body_style))
    elements.append(Spacer(1, 8))

    for i in range(num_target_genes):
        gene_label = f"Gene {i+1}"
        gene_rows_pdf = [r for r in stat_rows if r["__gene__"] == gene_label]
        if not gene_rows_pdf:
            continue
        try:
            fig_d, ax_d = plt.subplots(figsize=(6, 3.2))
            all_vals, all_labels = [list(gene_rows_pdf[0]["__ratio_ctrl__"])], ["Control"]
            for gr in gene_rows_pdf:
                all_vals.append(list(gr["__ratio_smp__"]))
                all_labels.append(gr["__group__"])
            palette = ['#00796b', '#d84315', '#455a64', '#f9a825', '#5e35b1']
            for k, (vals, lbl) in enumerate(zip(all_vals, all_labels)):
                col = palette[k % len(palette)]
                jitter = np.random.uniform(-0.08, 0.08, len(vals))
                ax_d.scatter([k + 1 + j for j in jitter], vals, color=col, alpha=0.75, s=28, zorder=3)
                ax_d.hlines(np.mean(vals), k + 0.75, k + 1.25, colors='black', linewidths=2, zorder=4)
            ax_d.set_xticks(range(1, len(all_labels) + 1))
            ax_d.set_xticklabels(all_labels, fontsize=8)
            ax_d.set_ylabel('Normalized Ratio', fontsize=9)
            ax_d.set_title(f'{gene_label} — Ratio Distribution', fontsize=10, fontweight='bold')
            ax_d.spines['top'].set_visible(False)
            ax_d.spines['right'].set_visible(False)
            plt.tight_layout()
            ib3 = BytesIO()
            plt.savefig(ib3, format='png', dpi=150, bbox_inches='tight')
            plt.close()
            ib3.seek(0)
            elements.append(RLImage(ib3, width=420, height=210))
            elements.append(Paragraph(safe_str(f"Figure. Ratio distribution for {gene_label}."), caption_style))
            elements.append(Spacer(1, 10))
        except Exception:
            pass
    elements.append(PageBreak())

    # ── SECTION 6: INTERPRETATION ─────────────────────────────────────────────
    elements.append(Paragraph(safe_str("6. How to Interpret Your Results"), h1_style))
    elements.append(hr())
    fc_hdr = ["Fold Change", "Interpretation", "Biological Significance"]
    fc_rows = [
        [">2.0", "Strong upregulation / gain", "Consider biologically relevant"],
        ["1.5-2.0", "Moderate upregulation / gain", "May be relevant; verify"],
        ["1.0-1.5", "Weak upregulation", "Likely not significant alone"],
        ["1.0", "No change", "No differential expression / CNV"],
        ["0.67-1.0", "Weak downregulation", "Likely not significant alone"],
        ["0.5-0.67", "Moderate downregulation / loss", "May be relevant; verify"],
        ["<0.5", "Strong downregulation / loss", "Consider biologically relevant"],
    ]
    elements.append(make_table([fc_hdr] + fc_rows, col_widths=[(letter[0] - 100) / 3] * 3))
    elements.append(Spacer(1, 8))
    elements.append(Paragraph(safe_str(
        "Note: statistical significance (p < 0.05) and biological significance (fold change "
        "magnitude, and the width of the Poisson 95% CI relative to the observed effect) "
        "should be considered together."
    ), warn_style))
    elements.append(PageBreak())

    # ── SECTION 7: REFERENCES ─────────────────────────────────────────────────
    elements.append(Paragraph(safe_str("7. References"), h1_style))
    elements.append(hr())
    refs = [
        "Yalçınkaya B (2026). GeneQuantify: a web-based tool for qPCR gene expression and copy number variation analysis. Molecular and Cellular Biochemistry. https://doi.org/10.1007/s11010-026-05621-y",
        "Huggett JF et al. (2013). The digital MIQE guidelines. Clinical Chemistry, 59(6), 892-902.",
        "Huggett JF (2020). The digital MIQE guidelines update. Clinical Chemistry, 66(8), 1012-1029.",
        "Vandesompele J et al. (2002). Genome Biology, 3(7). (geNorm normalization concept, adapted here to λ values)",
        "Grubbs FE (1969). Technometrics, 11(1), 1-21. (Outlier detection)",
        "Tukey JW (1977). Exploratory Data Analysis. Addison-Wesley. (IQR outlier method)",
        "Benjamini Y & Hochberg Y (1995). J Royal Stat Soc B, 57(1), 289-300. (FDR correction)",
        "Dube S, Qin J, Ramakrishnan R (2008). PLoS ONE, 3(8), e2876. (Poisson statistics for digital PCR)",
    ]
    for ref in refs:
        elements.append(Paragraph(safe_str(f"• {ref}"), small_style))
        elements.append(Spacer(1, 3))

    elements.append(Spacer(1, 16))
    elements.append(hr())
    elements.append(Paragraph(safe_str(
        f"AbsoluteGene — For research and educational use only. Not validated for clinical "
        f"diagnostic purposes. | Generated: {now} | Contact: mailtoburhanettin@gmail.com"
    ), small_style))

    doc.build(elements)
    buffer.seek(0)
    return buffer


with tab_report:
    st.markdown(f"### 📄 {_t['pdf_report']}")
    st.markdown("---")
    if not input_values_table:
        st.info(_t['error_no_data'])
    else:
        st.success('✅ ' + _t['pdf_ready'].format(n=len(input_values_table)))
        if st.button(f"📥 {_t['generate_pdf']}", key="pdf_btn"):
            pdf_buffer = create_pdf(data, stats_data, pd.DataFrame(input_values_table), language_code,
                                     multigroup_results=multigroup_results)
            st.download_button(
                label=f"⬇️ {_t['pdf_report']}", data=pdf_buffer,
                file_name="absolutegene_report.pdf", mime="application/pdf", key="pdf_dl"
            )

st.markdown(
    f"<h4 style='font-size: 12px; font-family: Arial, sans-serif; color: #555;'>"
    f"<a href='mailto:mailtoburhanettin@gmail.com' style='color: #555; text-decoration: none;'>"
    f"{_t['subtitle']}</a></h4>", unsafe_allow_html=True
)
