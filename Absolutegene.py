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
import json

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
.block-container { padding-top: 3rem !important; padding-bottom: 1rem !important; }
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
        "sidebar_example_title": "📋 Örnek Veri Yükle",
        "sidebar_example_select": "Senaryo seçin",
        "sidebar_example_load_btn": "▶ Senaryoyu Yükle",
        "sidebar_example_loaded": "✅ {s} yüklendi! Veri Girişi sekmesine geçin.",
        "ntc_expander": "🧫 NTC (Negatif Kontrol) / LOD-LOQ — opsiyonel",
        "ntc_description": "Tespit sınırını (LOD) ve ölçüm sınırını (LOQ) hesaplamak için No-Template Control (şablon içermeyen kontrol) replikatlarınızı girin. Boş bırakılırsa bu adım atlanır.",
        "ntc_positive_label": "NTC Pozitif Partisyon Sayısı",
        "ntc_total_label": "NTC Toplam Kabul Edilen Partisyon Sayısı",
        "ntc_calc_btn": "📊 LOD/LOQ Hesapla",
        "lod_result_title": "Tespit Sınırı (LOD) ve Ölçüm Sınırı (LOQ) — {gene}",
        "lod_label": "LOD (kopya/µL)",
        "loq_label": "LOQ (kopya/µL)",
        "ntc_zero_note": "ℹ️ NTC replikatlarında hiç pozitif partisyon tespit edilmedi (n={n} partisyon havuzlandı). LOD, 'üçler kuralı' (rule of three) ile üst güven sınırından hesaplandı: 3/n.",
        "ntc_contamination_warning": "⚠️ NTC replikatlarında arka plan sinyali tespit edildi (λ={lam:.5f}, {pos}/{tot} pozitif partisyon). Bu, reaktif kontaminasyonu veya çapraz kontaminasyon olabilir — LOD, bu arka plan seviyesi + güvenlik payı olarak hesaplandı.",
        "below_lod_flag": "❌ LOD altı",
        "between_lod_loq_flag": "⚠️ LOD-LOQ arası",
        "above_loq_flag": "✅ LOQ üstü",
        "lod_qc_col": "LOD/LOQ Durumu",
        "loq_heuristic_note": "Not: LOQ, yaygın kullanılan bir kural olarak LOD'un 3 katı şeklinde hesaplanmıştır. Daha kesin bir LOQ için tekrarlanabilirlik (CV%) temelli deneysel doğrulama önerilir.",
        "csv_import_expander": "📂 Cihaz CSV Dosyası İçe Aktar (QuantaSoft / QX Manager / QIAcuity)",
        "csv_import_description": "Bio-Rad QuantaSoft/QX Manager veya QIAGEN QIAcuity yazılımından dışa aktardığınız CSV dosyasını yükleyin. Sütunlar otomatik tanınmaya çalışılır; gerekirse manuel eşleştirme yapabilirsiniz.",
        "csv_uploader": "CSV dosyası seçin",
        "csv_parse_error": "❌ Dosya ayrıştırma hatası: {err}",
        "csv_col_mapping_title": "**Sütun Eşleştirme**",
        "csv_col_sample": "Örnek Adı Sütunu",
        "csv_col_target": "Hedef/Assay Adı Sütunu",
        "csv_col_positives": "Pozitif Partisyon Sütunu",
        "csv_col_total": "Toplam Partisyon Sütunu",
        "csv_col_auto_detected": "✅ Otomatik tespit edildi",
        "csv_col_manual_needed": "⚠️ Otomatik tespit edilemedi — lütfen manuel seçin",
        "csv_preview_title": "Ayrıştırılan verileri önizle",
        "csv_assay_assignment_title": "**Assay (Hedef/Referans Gen) Ataması**",
        "csv_target_assays_label": "Hedef gen assay(ler)i",
        "csv_ref_assays_label": "Referans gen assay(ler)i",
        "csv_group_assignment_title": "**Grup Ataması**",
        "csv_ctrl_label": "Kontrol örnek adı (virgülle ayrılmış alt dizeler)",
        "csv_n_patient_groups": "Hasta grubu sayısı",
        "csv_patient_label": "Hasta grubu {i} örnek adı (adları)",
        "csv_apply_btn": "✅ İçe Aktarımı Veri Girişine Uygula",
        "csv_apply_success": "✅ {n} değer Veri Girişi sekmesine yüklendi! Kontrol edip ayarlayabilirsiniz.",
        "csv_apply_warning": "⚠️ Hiçbir değer eşleştirilemedi. Kontrol/hasta etiketlerinin örnek adlarıyla uyuştuğundan emin olun.",
        "download_example_csv": "📄 Örnek CSV Dosyasını İndir",
        "dilution_expander": "🧪 Optimal Seyreltme / Dinamik Aralık Hesaplayıcısı",
        "dilution_description": "Poisson tabanlı kantifikasyonun en hassas olduğu aralık partisyon başına ~1.0–2.0 kopyadır (λ). Bu hesaplayıcı, mevcut ölçümünüze göre önerilen seyreltme faktörünü bulur.",
        "dilution_mode_label": "Giriş yöntemi",
        "dilution_mode_counts": "Gözlenen partisyon sayılarından",
        "dilution_mode_lambda": "Bilinen λ değerinden",
        "dilution_positive_label": "Pozitif Partisyon Sayısı",
        "dilution_total_label": "Toplam Partisyon Sayısı",
        "dilution_lambda_label": "Mevcut λ (kopya/partisyon)",
        "dilution_target_label": "Hedef λ (kopya/partisyon)",
        "dilution_calc_btn": "📊 Seyreltme Öner",
        "dilution_result_optimal": "✅ **Optimal aralıkta** (λ={lam:.3f}). Ek seyreltme gerekmez.",
        "dilution_result_too_low": "⚠️ **Çok seyreltik** (λ={lam:.3f}). Bir sonraki hazırlığınızda seyreltme faktörünüzü **mevcut faktörün {factor:.3f} katı** yapın (yani daha az seyreltin).",
        "dilution_result_too_high": "⚠️ **Optimal aralığın üzerinde** (λ={lam:.3f}). Bir sonraki hazırlığınızda seyreltme faktörünüzü **mevcut faktörün {factor:.3f} katı** yapın (yani daha fazla seyreltin).",
        "dilution_result_saturated": "❌ **Doygunluğa yakın/riskli** (λ={lam:.3f}). Bir sonraki hazırlığınızda seyreltme faktörünüzü **mevcut faktörün {factor:.3f} katı** yapın.",
        "dilution_factor_note": "Bu hesaplayıcı, orijinal seyreltme faktörünüzü bilmediği için mutlak bir seyreltme oranı değil, **göreli bir çarpan** önerir: bir sonraki hazırlığınızda kullandığınız seyreltme faktörünü bu çarpanla çarpın. Örnek: şu an 1:50 seyreltiyorsanız ve önerilen çarpan 2.0 ise, bir sonraki denemede 1:100 seyreltin.",
        "qc_panel_title": "🩺 Veri Kalitesi Özeti",
        "qc_panel_total": "Toplam Replikat",
        "qc_panel_qc_fail": "Düşük Partisyon (QC)",
        "qc_panel_saturated": "Doygun",
        "qc_panel_outlier": "Aykırı Değer (Dışlanan)",
        "qc_panel_below_lod": "LOD Altı Sonuç",
        "qc_panel_high_cv": "Yüksek CV (>%25)",
        "qc_panel_verdict_good": "✅ **İyi** — veri seti sağlıklı görünüyor, önemli bir kalite sorunu tespit edilmedi.",
        "qc_panel_verdict_caution": "⚠️ **Dikkat** — bazı replikatlarda kalite bayrakları var. Sonuçları yorumlarken bu noktaları göz önünde bulundurun.",
        "qc_panel_verdict_poor": "❌ **Gözden geçirin** — veri setinde önemli oranda kalite sorunu bayrağı var (QC hatası, doygunluk veya yüksek değişkenlik). Rapor öncesi verileri kontrol etmenizi öneririz.",
        "qc_panel_no_data": "Henüz veri girilmedi.",
        "project_expander": "💾 Proje Kaydet / Yükle",
        "project_description": "Tüm veri girişinizi (çalışma tasarımı, kontrol/hasta grubu verileri, NTC) bir JSON dosyası olarak kaydedip daha sonra kaldığınız yerden devam edebilirsiniz.",
        "project_export_btn": "📥 Projeyi Dışa Aktar (JSON)",
        "project_import_uploader": "Proje dosyası yükleyin (.json)",
        "project_import_success": "✅ Proje yüklendi ({n} alan geri yüklendi)! Veri Girişi sekmesine geçin.",
        "project_import_error": "❌ Proje dosyası okunamadı: {err}",
        "batch_pdf_btn": "📥 Tarama Raporunu PDF Olarak Hazırla",
        "batch_pdf_report": "Toplu Numune Tarama Raporu",
        "batch_pdf_description": "Bu rapor, her numune için Poisson istatistiğine dayalı λ, normalize oran, %95 güven aralığı ve beklenen orana göre sınıflandırmayı içerir. Sınıflandırma, örneğin oran güven aralığının beklenen (değişim-yok) değeri içerip içermediğine göre belirlenir.",
        "vaf_pdf_btn": "📥 VAF Raporunu PDF Olarak Hazırla",
        "vaf_pdf_report": "VAF / Mutasyon Fraksiyonu Raporu",
        "vaf_pdf_description": "Bu rapor, her örnek için mutant ve wild-type assay λ değerlerini, Fraksiyonel Bolluk (%FA) ve delta-method %95 güven aralığını içerir.",
        "tab_clinical": "Klinik Araçlar",
        "tab_crm": "CRM Üretimi",
        "clinical_title": "🩺 Klinik Doğrulama Araçları",
        "clinical_description": "Ölçüm belirsizliği, referans değişim değeri, hassasiyet çalışması ve yöntem karşılaştırması gibi klinik/tanısal doğrulama için gerekli araçlar. Bu araçlar araştırma ve yöntem doğrulama amaçlıdır; klinik tanı kararı için doğrulanmamıştır.",
        "clinical_mode_label": "Araç seçin",
        "clinical_mode_mu": "📐 Ölçüm Belirsizliği (MU) Bütçesi",
        "clinical_mode_rcv": "📈 Referans Değişim Değeri (RCV)",
        "clinical_mode_precision": "🔁 Hassasiyet Çalışması",
        "clinical_mode_comparison": "⚖️ Yöntem Karşılaştırma",
        # MU Budget
        "mu_description": "GUM (Ölçüm Belirsizliği İfadesi Kılavuzu) yaklaşımıyla, bağımsız belirsizlik kaynakları karesel olarak birleştirilir ve bir kapsam faktörü (k) ile çarpılarak genişletilmiş belirsizlik elde edilir.",
        "mu_poisson_source": "Poisson Sayım İstatistiği Belirsizliği (%)",
        "mu_poisson_help": "Otomatik hesaplanabilir (replikat CV%'den) veya manuel girilebilir.",
        "mu_poisson_auto": "Sonuçlardan otomatik al",
        "mu_poisson_manual": "Manuel gir",
        "mu_pipetting_label": "Pipetleme / Seyreltme Belirsizliği (%)",
        "mu_pipetting_help": "Tipik olarak iyi pipetleme tekniğinde %1-2 CV. Kendi doğrulama verinize göre girin.",
        "mu_precision_label": "Run-arası Hassasiyet Belirsizliği (%, opsiyonel)",
        "mu_precision_help": "Ayrı bir hassasiyet çalışmasından elde edilen günler-arası CV%. Bilinmiyorsa 0 bırakın.",
        "mu_gene_select": "Sonuç seçin (Gen / Grup)",
        "mu_k_label": "Kapsam Faktörü (k)",
        "mu_calc_btn": "📊 Belirsizlik Bütçesini Hesapla",
        "mu_result_title": "Sonuç: {value} ± {U:.2f}% (k={k}, ~%95 güven düzeyi)",
        "mu_budget_table_title": "**Belirsizlik Bütçesi Tablosu**",
        "mu_col_source": "Kaynak",
        "mu_col_contribution": "Katkı (%, göreli standart belirsizlik)",
        "mu_col_combined": "Birleşik Standart Belirsizlik (u_c)",
        "mu_col_expanded": "Genişletilmiş Belirsizlik (U)",
        # RCV
        "rcv_description": "İki seri sonuç arasındaki farkın, analitik ve biyolojik değişkenlikten beklenenin ötesinde istatistiksel olarak anlamlı olup olmadığını belirler (Fraser & Harris, Crit Rev Clin Lab Sci 1989).",
        "rcv_result1_label": "Sonuç 1 (önceki, örn. %VAF veya kopya/µL)",
        "rcv_result2_label": "Sonuç 2 (sonraki)",
        "rcv_cv_analytical_label": "Analitik CV (%)",
        "rcv_cv_analytical_help": "Bu ölçümün tekrarlanabilirlik CV%'si (replikat CV veya MU bütçesinden).",
        "rcv_cv_biological_label": "Biyolojik Değişkenlik CV (%)",
        "rcv_cv_biological_help": "Analitin bilinen birey-içi biyolojik değişkenliği (literatürden — analite özgüdür).",
        "rcv_z_label": "z-değeri (güven düzeyi)",
        "rcv_calc_btn": "📊 RCV Hesapla",
        "rcv_result_significant": "🔴 **Anlamlı değişim** — %{change:+.2f} değişim, RCV eşiği olan %{rcv:.2f}'yi aşıyor. Bu değişim ölçüm gürültüsüyle açıklanamaz.",
        "rcv_result_not_significant": "🟢 **Anlamlı değişim yok** — %{change:+.2f} değişim, RCV eşiği olan %{rcv:.2f} içinde. Analitik/biyolojik gürültü ile tutarlı.",
        # Precision study
        "precision_description": "Gün-içi (tekrarlanabilirlik) ve günler-arası (ara hassasiyet) varyans bileşenlerini iç içe (nested) ANOVA ile ayırır (CLSI EP05-A3 yaklaşımı).",
        "precision_n_days": "Gün/Çalışma Sayısı",
        "precision_day_label": "Gün {i} — Replikat Değerleri (her satıra bir değer)",
        "precision_calc_btn": "📊 Hassasiyeti Hesapla",
        "precision_repeatability_cv": "Tekrarlanabilirlik CV (%)",
        "precision_between_day_cv": "Günler-arası CV (%)",
        "precision_total_cv": "Toplam (Ara Hassasiyet) CV (%)",
        "precision_grand_mean": "Genel Ortalama",
        # Method comparison
        "comparison_description": "İki ölçüm yöntemi arasındaki uyumu Bland-Altman analizi ve Deming regresyonu ile değerlendirir.",
        "comparison_method1_label": "Yöntem 1 Değerleri (örn. dPCR — her satıra bir değer)",
        "comparison_method2_label": "Yöntem 2 Değerleri (örn. qPCR — eşleştirilmiş, aynı sırada)",
        "comparison_variance_ratio_label": "Varyans Oranı (λ, Deming için)",
        "comparison_calc_btn": "📊 Karşılaştırmayı Hesapla",
        "comparison_bias_label": "Ortalama Fark (Bias)",
        "comparison_loa_label": "Uyum Sınırları (%95)",
        "comparison_deming_label": "Deming Regresyonu",
        "comparison_ba_chart_title": "Bland-Altman Grafiği",
        "comparison_deming_chart_title": "Deming Regresyon Grafiği",
        # CRM Production
        "crm_title": "🏭 Sertifikalı Referans Malzeme (CRM) Üretim Araçları",
        "crm_description": "ISO Guide 35 / ISO 17034 genel ilkelerine dayalı homojenlik testi, stabilite testi, atanmış değer & belirsizlik bütçesi ve sertifika (CoA) oluşturucu. Bu araçlar araştırma amaçlıdır; resmi bir CRM üretim/onay süreci yerine geçmez.",
        "crm_mode_label": "Araç seçin",
        "crm_mode_homogeneity": "🧪 Homojenlik Testi",
        "crm_mode_stability": "⏳ Stabilite Testi",
        "crm_mode_uncertainty": "📐 Atanmış Değer & Belirsizlik Bütçesi",
        "crm_mode_coa": "📜 Sertifika (CoA) Oluşturucu",
        # Homogeneity
        "homog_description": "Bir üretim partisindeki çok sayıda ünite/viyal arasında tek yönlü ANOVA ile homojenlik testi (Linsinger ve ark. 2001). F ≤ F_kritik ise partinin homojen olduğu kabul edilir.",
        "homog_n_units": "Ünite/Viyal Sayısı",
        "homog_unit_label": "Ünite {i} — Replikat Değerleri (her satıra bir değer)",
        "homog_calc_btn": "📊 Homojenliği Test Et",
        "homog_result_homogeneous": "✅ **Homojen** (F={F:.3f} ≤ F_kritik={fcrit:.3f}, p={p:.4f}). Partide istatistiksel olarak anlamlı bir üniteler-arası fark tespit edilmedi.",
        "homog_result_inhomogeneous": "❌ **Homojen Değil** (F={F:.3f} > F_kritik={fcrit:.3f}, p={p:.4f}). Üniteler arasında istatistiksel olarak anlamlı fark var — parti gözden geçirilmeli.",
        "homog_ubb_label": "Homojenlik Belirsizlik Katkısı (u_bb)",
        "homog_grand_mean_label": "Genel Ortalama",
        # Stability
        "stab_description": "Zaman içindeki ölçüm değerlerine doğrusal regresyon uygulayarak anlamlı bir bozunma trendi olup olmadığını test eder (ISO Guide 35 klasik stabilite yaklaşımı).",
        "stab_time_label": "Zaman Noktaları (örn. gün/ay — her satıra bir değer)",
        "stab_value_label": "Ölçülen Değerler (zaman noktalarıyla eşleştirilmiş, aynı sırada)",
        "stab_duration_label": "İncelenen Süre (raf ömrü, zaman noktalarıyla aynı birimde)",
        "stab_calc_btn": "📊 Stabiliteyi Test Et",
        "stab_result_stable": "✅ **Kararlı** (eğim p={p:.4f} ≥ 0.05). Zaman içinde istatistiksel olarak anlamlı bir trend tespit edilmedi.",
        "stab_result_unstable": "⚠️ **Anlamlı Trend Tespit Edildi** (eğim p={p:.4f} < 0.05). Materyal zaman içinde değişim gösteriyor olabilir — raf ömrü iddiasını gözden geçirin.",
        "stab_ustab_label": "Stabilite Belirsizlik Katkısı (u_stab)",
        "stab_slope_label": "Eğim (birim/zaman)",
        "stab_chart_title": "Stabilite — Zamana Karşı Değer",
        # Uncertainty budget
        "unc_description": "Karakterizasyon, homojenlik ve stabilite belirsizlik bileşenlerini GUM yaklaşımıyla (karesel toplam) birleştirerek nihai genişletilmiş belirsizliği hesaplar. Homojenlik/Stabilite sekmelerinde hesaplama yaptıysanız değerler otomatik doldurulur.",
        "unc_assigned_value_label": "Atanmış Değer",
        "unc_u_char_label": "Karakterizasyon Belirsizliği (u_char, mutlak birim)",
        "unc_u_char_help": "Atanmış değeri belirlemek için kullanılan ölçümlerin standart belirsizliği (örn. ortalamanın standart hatası).",
        "unc_u_bb_label": "Homojenlik Belirsizliği (u_bb)",
        "unc_u_stab_label": "Stabilite Belirsizliği (u_stab)",
        "unc_k_label": "Kapsam Faktörü (k)",
        "unc_use_cached": "Homojenlik/Stabilite sekmelerinden otomatik doldur",
        "unc_calc_btn": "📊 Belirsizlik Bütçesini Hesapla",
        "unc_result_title": "Atanmış Değer = {value:.5f} ± {U:.5f} (k={k}, %{urel:.2f} göreli)",
        "unc_no_cache": "ℹ️ Önceden hesaplanmış Homojenlik/Stabilite sonucu bulunamadı — manuel giriş kullanılacak.",
        # CoA
        "coa_description": "Atanmış değer, genişletilmiş belirsizlik ve izlenebilirlik bilgilerini içeren resmi tarzda bir Analiz Sertifikası (CoA) PDF'i oluşturur.",
        "coa_material_name": "Materyal Adı",
        "coa_lot_number": "Parti/Lot Numarası",
        "coa_producer": "Üretici/Laboratuvar",
        "coa_assigned_value_label": "Atanmış Değer",
        "coa_unit_label": "Birim",
        "coa_expanded_unc_label": "Genişletilmiş Belirsizlik (U)",
        "coa_k_label": "Kapsam Faktörü (k)",
        "coa_traceability_label": "İzlenebilirlik İfadesi",
        "coa_traceability_default": "Bu değerin izlenebilirliği, dijital PCR ile Poisson istatistiğine dayalı mutlak kantifikasyon yoluyla SI birimine (mol) dayanmaktadır.",
        "coa_validity_label": "Geçerlilik Tarihi / Raf Ömrü",
        "coa_generate_btn": "📥 Sertifikayı Oluştur (PDF)",
        "coa_download_btn": "⬇️ Analiz Sertifikası (CoA)",
        "tab_batch": "Toplu Tarama",
        "batch_title": "🔬 Toplu Numune Tarama (CNV Taraması)",
        "batch_description": "Çok sayıda numuneyi (örn. bir kohort) tek bir referans/beklenen orana karşı hızlıca taramak için tasarlanmıştır. Her numune için ayrı replikat grupları oluşturmak yerine, CSV dosyanızı yükleyip her örnek için tek bir Poisson tabanlı güven aralığı hesaplanır.",
        "batch_uploader": "Toplu tarama CSV dosyası yükleyin",
        "batch_assay_title": "**Assay Ataması**",
        "batch_target_label": "Hedef gen assay'i",
        "batch_ref_label": "Referans gen assay(ler)i (birden fazla seçilirse partisyonlar havuzlanır)",
        "batch_settings_title": "**Tarama Ayarları**",
        "batch_expected_ratio": "Beklenen Oran (değişim yok referansı)",
        "batch_ploidy": "Referans Ploidi",
        "batch_run_btn": "▶ Taramayı Çalıştır",
        "batch_results_title": "📋 Tarama Sonuçları",
        "batch_col_sample": "Örnek",
        "batch_col_lambda_t": "λ Hedef",
        "batch_col_lambda_r": "λ Referans",
        "batch_col_ratio": "Oran",
        "batch_col_ci": "%95 GA (Oran)",
        "batch_col_cn": "Kopya Sayısı (CN)",
        "batch_col_conc": "Konsantrasyon (kopya/µL)",
        "batch_col_class": "Sınıflandırma",
        "batch_class_gain": "📈 Kazanım",
        "batch_class_loss": "📉 Kayıp",
        "batch_class_normal": "➖ Normal",
        "batch_chart_title": "Numuneler Arası Kopya Sayısı Taraması",
        "batch_download_btn": "📥 Tarama Sonuçlarını İndir (CSV)",
        "batch_no_data": "Sonuç yok — önce dosya yükleyip taramayı çalıştırın.",
        "batch_n_flagged": "{n} / {total} numune beklenen orandan istatistiksel olarak anlamlı şekilde farklı (%95 GA değişim-yok değerini içermiyor).",
        "tab_vaf": "VAF Hesaplayıcı",
        "vaf_title": "🧬 VAF / Mutasyon Fraksiyonu Hesaplayıcı",
        "vaf_description": "Likit biyopsi / ctDNA izleme için tasarlanmıştır. Mutant ve Wild-type (yabanıl tip) assay'lerini içeren bir CSV yükleyin; her örnek için Fraksiyonel Bolluk (FA%, varyant allel fraksiyonu) ve %95 güven aralığı hesaplanır.",
        "vaf_uploader": "VAF CSV dosyası yükleyin",
        "vaf_assay_title": "**Assay Ataması**",
        "vaf_mutant_label": "Mutant assay",
        "vaf_wt_label": "Wild-type (yabanıl tip) assay",
        "vaf_run_btn": "▶ VAF Hesapla",
        "vaf_results_title": "📋 VAF Sonuçları",
        "vaf_col_sample": "Örnek",
        "vaf_col_lambda_mut": "λ Mutant",
        "vaf_col_lambda_wt": "λ Wild-type",
        "vaf_col_fa": "FA (%)",
        "vaf_col_ci": "%95 GA (FA%)",
        "vaf_col_conc_mut": "Mutant Konsantrasyonu (kopya/µL)",
        "vaf_col_detected": "Tespit Durumu",
        "vaf_detected_yes": "✅ Tespit Edildi",
        "vaf_detected_no": "❌ Tespit Edilemedi",
        "vaf_chart_title": "Örnekler Arası Fraksiyonel Bolluk (VAF%)",
        "vaf_download_btn": "📥 VAF Sonuçlarını İndir (CSV)",
        "vaf_no_data": "Sonuç yok — önce dosya yükleyip hesaplamayı çalıştırın.",
        "vaf_method_note": "ℹ️ FA = λ(mutant) / (λ(mutant) + λ(wild-type)). %95 GA delta-method ile hesaplanmıştır (Hindson et al. 2013, Anal Chem). Bu hesaplama, mutant ve wild-type assay'lerinin aynı partisyon setinde (multipleks) veya eşdeğer koşullarda ölçüldüğünü varsayar.",
        "vaf_n_detected": "{n} / {total} örnekte mutant alel tespit edildi (≥1 pozitif partisyon).",
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
        "sidebar_example_title": "📋 Load Example Data",
        "sidebar_example_select": "Select scenario",
        "sidebar_example_load_btn": "▶ Load Scenario",
        "sidebar_example_loaded": "✅ {s} loaded! Switch to the Data Entry tab.",
        "ntc_expander": "🧫 NTC (No-Template Control) / LOD-LOQ — optional",
        "ntc_description": "Enter your No-Template Control replicates to calculate the limit of detection (LOD) and limit of quantification (LOQ). Leave blank to skip this step.",
        "ntc_positive_label": "NTC Positive Partition Count",
        "ntc_total_label": "NTC Total Accepted Partition Count",
        "ntc_calc_btn": "📊 Calculate LOD/LOQ",
        "lod_result_title": "Limit of Detection (LOD) and Limit of Quantification (LOQ) — {gene}",
        "lod_label": "LOD (copies/µL)",
        "loq_label": "LOQ (copies/µL)",
        "ntc_zero_note": "ℹ️ No positive partitions detected across NTC replicates (n={n} partitions pooled). LOD was calculated using the rule-of-three upper confidence bound: 3/n.",
        "ntc_contamination_warning": "⚠️ Background signal detected in NTC replicates (λ={lam:.5f}, {pos}/{tot} positive partitions). This may indicate reagent or cross-contamination — LOD was calculated as this background level plus a safety margin.",
        "below_lod_flag": "❌ Below LOD",
        "between_lod_loq_flag": "⚠️ Between LOD-LOQ",
        "above_loq_flag": "✅ Above LOQ",
        "lod_qc_col": "LOD/LOQ Status",
        "loq_heuristic_note": "Note: LOQ is calculated as 3× LOD, a commonly used heuristic. For a more rigorous LOQ, empirical validation based on reproducibility (CV%) is recommended.",
        "csv_import_expander": "📂 Import Instrument CSV (QuantaSoft / QX Manager / QIAcuity)",
        "csv_import_description": "Upload a CSV export from Bio-Rad QuantaSoft/QX Manager or QIAGEN QIAcuity software. Columns are auto-detected where possible; you can manually map any that aren't recognized.",
        "csv_uploader": "Choose CSV file",
        "csv_parse_error": "❌ File parse error: {err}",
        "csv_col_mapping_title": "**Column Mapping**",
        "csv_col_sample": "Sample Name Column",
        "csv_col_target": "Target/Assay Name Column",
        "csv_col_positives": "Positive Partitions Column",
        "csv_col_total": "Total Partitions Column",
        "csv_col_auto_detected": "✅ Auto-detected",
        "csv_col_manual_needed": "⚠️ Not auto-detected — please select manually",
        "csv_preview_title": "Preview parsed data",
        "csv_assay_assignment_title": "**Assay (Target/Reference Gene) Assignment**",
        "csv_target_assays_label": "Target gene assay(s)",
        "csv_ref_assays_label": "Reference gene assay(s)",
        "csv_group_assignment_title": "**Group Assignment**",
        "csv_ctrl_label": "Control sample name(s) (comma-separated substrings)",
        "csv_n_patient_groups": "Number of patient groups",
        "csv_patient_label": "Patient group {i} sample name(s)",
        "csv_apply_btn": "✅ Apply Import to Data Entry",
        "csv_apply_success": "✅ {n} value(s) loaded into the Data Entry tab! Switch to review and adjust.",
        "csv_apply_warning": "⚠️ No values were mapped. Check that your control/patient labels match the sample names.",
        "download_example_csv": "📄 Download Example CSV",
        "dilution_expander": "🧪 Optimal Dilution / Dynamic Range Calculator",
        "dilution_description": "Poisson-based quantification is most precise in the range of ~1.0–2.0 copies per partition (λ). This calculator recommends a dilution factor based on your current measurement.",
        "dilution_mode_label": "Input method",
        "dilution_mode_counts": "From observed partition counts",
        "dilution_mode_lambda": "From known λ value",
        "dilution_positive_label": "Positive Partition Count",
        "dilution_total_label": "Total Partition Count",
        "dilution_lambda_label": "Current λ (copies/partition)",
        "dilution_target_label": "Target λ (copies/partition)",
        "dilution_calc_btn": "📊 Recommend Dilution",
        "dilution_result_optimal": "✅ **In the optimal range** (λ={lam:.3f}). No further dilution needed.",
        "dilution_result_too_low": "⚠️ **Too dilute** (λ={lam:.3f}). For your next prep, set your dilution factor to **{factor:.3f}× your current factor** (i.e., dilute less).",
        "dilution_result_too_high": "⚠️ **Above the optimal range** (λ={lam:.3f}). For your next prep, set your dilution factor to **{factor:.3f}× your current factor** (i.e., dilute more).",
        "dilution_result_saturated": "❌ **Near/at saturation risk** (λ={lam:.3f}). For your next prep, set your dilution factor to **{factor:.3f}× your current factor**.",
        "dilution_factor_note": "Since this calculator doesn't know your original dilution factor, it recommends a **relative multiplier**: multiply the dilution factor you used by this value for your next prep. Example: if you currently dilute 1:50 and the recommended factor is 2.0, try 1:100 next time.",
        "qc_panel_title": "🩺 Data Quality Summary",
        "qc_panel_total": "Total Replicates",
        "qc_panel_qc_fail": "Low Partition (QC)",
        "qc_panel_saturated": "Saturated",
        "qc_panel_outlier": "Outlier (Excluded)",
        "qc_panel_below_lod": "Below-LOD Results",
        "qc_panel_high_cv": "High CV (>25%)",
        "qc_panel_verdict_good": "✅ **Good** — the dataset looks healthy, no significant quality issues detected.",
        "qc_panel_verdict_caution": "⚠️ **Caution** — some replicates have quality flags. Consider these when interpreting results.",
        "qc_panel_verdict_poor": "❌ **Review recommended** — a significant proportion of quality flags (QC failures, saturation, or high variability) were detected. We recommend reviewing the data before reporting.",
        "qc_panel_no_data": "No data entered yet.",
        "project_expander": "💾 Save / Load Project",
        "project_description": "Save all your data entry (study design, control/patient group data, NTC) as a JSON file, then resume where you left off later.",
        "project_export_btn": "📥 Export Project (JSON)",
        "project_import_uploader": "Upload project file (.json)",
        "project_import_success": "✅ Project loaded ({n} fields restored)! Switch to the Data Entry tab.",
        "project_import_error": "❌ Could not read project file: {err}",
        "batch_pdf_btn": "📥 Prepare Screening Report (PDF)",
        "batch_pdf_report": "Batch Sample Screening Report",
        "batch_pdf_description": "This report includes, for each sample, the Poisson-derived λ, normalized ratio, 95% confidence interval, and classification relative to the expected ratio. Classification is based on whether the sample's ratio confidence interval includes the expected (no-change) value.",
        "vaf_pdf_btn": "📥 Prepare VAF Report (PDF)",
        "vaf_pdf_report": "VAF / Mutation Fraction Report",
        "vaf_pdf_description": "This report includes, for each sample, the mutant and wild-type assay λ values, Fractional Abundance (FA%), and the delta-method 95% confidence interval.",
        "tab_clinical": "Clinical Tools",
        "tab_crm": "CRM Production",
        "clinical_title": "🩺 Clinical Validation Tools",
        "clinical_description": "Tools for clinical/diagnostic validation: measurement uncertainty, reference change value, precision study, and method comparison. These tools are for research and method-validation purposes; not validated for clinical diagnostic decision-making.",
        "clinical_mode_label": "Select tool",
        "clinical_mode_mu": "📐 Measurement Uncertainty (MU) Budget",
        "clinical_mode_rcv": "📈 Reference Change Value (RCV)",
        "clinical_mode_precision": "🔁 Precision Study",
        "clinical_mode_comparison": "⚖️ Method Comparison",
        "mu_description": "Using the GUM (Guide to the Expression of Uncertainty in Measurement) approach, independent uncertainty sources are combined in quadrature and multiplied by a coverage factor (k) to obtain the expanded uncertainty.",
        "mu_poisson_source": "Poisson Counting Statistics Uncertainty (%)",
        "mu_poisson_help": "Can be auto-calculated (from replicate CV%) or entered manually.",
        "mu_poisson_auto": "Auto-fetch from results",
        "mu_poisson_manual": "Enter manually",
        "mu_pipetting_label": "Pipetting / Dilution Uncertainty (%)",
        "mu_pipetting_help": "Typically ~1-2% CV for good pipetting technique. Enter based on your own validation data.",
        "mu_precision_label": "Inter-run Precision Uncertainty (%, optional)",
        "mu_precision_help": "Between-day CV% from a separate precision study. Leave at 0 if unknown.",
        "mu_gene_select": "Select result (Gene / Group)",
        "mu_k_label": "Coverage Factor (k)",
        "mu_calc_btn": "📊 Calculate Uncertainty Budget",
        "mu_result_title": "Result: {value} ± {U:.2f}% (k={k}, ~95% confidence level)",
        "mu_budget_table_title": "**Uncertainty Budget Table**",
        "mu_col_source": "Source",
        "mu_col_contribution": "Contribution (%, relative standard uncertainty)",
        "mu_col_combined": "Combined Standard Uncertainty (u_c)",
        "mu_col_expanded": "Expanded Uncertainty (U)",
        "rcv_description": "Determines whether the difference between two serial results is statistically significant beyond what would be expected from analytical and biological variability alone (Fraser & Harris, Crit Rev Clin Lab Sci 1989).",
        "rcv_result1_label": "Result 1 (earlier, e.g. %VAF or copies/µL)",
        "rcv_result2_label": "Result 2 (later)",
        "rcv_cv_analytical_label": "Analytical CV (%)",
        "rcv_cv_analytical_help": "The repeatability CV% of this measurement (replicate CV or from the MU budget).",
        "rcv_cv_biological_label": "Biological Variation CV (%)",
        "rcv_cv_biological_help": "The known within-subject biological variation of the analyte (from literature — analyte-specific).",
        "rcv_z_label": "z-value (confidence level)",
        "rcv_calc_btn": "📊 Calculate RCV",
        "rcv_result_significant": "🔴 **Significant change** — {change:+.2f}% change exceeds the RCV threshold of {rcv:.2f}%. This change cannot be explained by measurement noise alone.",
        "rcv_result_not_significant": "🟢 **No significant change** — {change:+.2f}% change is within the RCV threshold of {rcv:.2f}%. Consistent with analytical/biological noise.",
        "precision_description": "Separates within-day (repeatability) and between-day (intermediate precision) variance components using nested ANOVA (CLSI EP05-A3 approach).",
        "precision_n_days": "Number of Days/Runs",
        "precision_day_label": "Day {i} — Replicate Values (one per line)",
        "precision_calc_btn": "📊 Calculate Precision",
        "precision_repeatability_cv": "Repeatability CV (%)",
        "precision_between_day_cv": "Between-Day CV (%)",
        "precision_total_cv": "Total (Intermediate Precision) CV (%)",
        "precision_grand_mean": "Grand Mean",
        "comparison_description": "Assesses agreement between two measurement methods using Bland-Altman analysis and Deming regression.",
        "comparison_method1_label": "Method 1 Values (e.g. dPCR — one per line)",
        "comparison_method2_label": "Method 2 Values (e.g. qPCR — paired, same order)",
        "comparison_variance_ratio_label": "Variance Ratio (λ, for Deming)",
        "comparison_calc_btn": "📊 Calculate Comparison",
        "comparison_bias_label": "Mean Difference (Bias)",
        "comparison_loa_label": "Limits of Agreement (95%)",
        "comparison_deming_label": "Deming Regression",
        "comparison_ba_chart_title": "Bland-Altman Plot",
        "comparison_deming_chart_title": "Deming Regression Plot",
        # CRM Production
        "crm_title": "🏭 Certified Reference Material (CRM) Production Tools",
        "crm_description": "Homogeneity testing, stability testing, assigned value & uncertainty budget, and certificate (CoA) generation, based on the general principles of ISO Guide 35 / ISO 17034. These tools are for research purposes and do not replace a formal CRM production/certification process.",
        "crm_mode_label": "Select tool",
        "crm_mode_homogeneity": "🧪 Homogeneity Testing",
        "crm_mode_stability": "⏳ Stability Testing",
        "crm_mode_uncertainty": "📐 Assigned Value & Uncertainty Budget",
        "crm_mode_coa": "📜 Certificate (CoA) Generator",
        "homog_description": "Tests homogeneity across multiple units/vials from a production batch using one-way ANOVA (Linsinger et al. 2001). The batch is considered homogeneous if F ≤ F_critical.",
        "homog_n_units": "Number of Units/Vials",
        "homog_unit_label": "Unit {i} — Replicate Values (one per line)",
        "homog_calc_btn": "📊 Test Homogeneity",
        "homog_result_homogeneous": "✅ **Homogeneous** (F={F:.3f} ≤ F_critical={fcrit:.3f}, p={p:.4f}). No statistically significant between-unit difference was detected.",
        "homog_result_inhomogeneous": "❌ **Not Homogeneous** (F={F:.3f} > F_critical={fcrit:.3f}, p={p:.4f}). A statistically significant between-unit difference was found — the batch should be reviewed.",
        "homog_ubb_label": "Homogeneity Uncertainty Contribution (u_bb)",
        "homog_grand_mean_label": "Grand Mean",
        "stab_description": "Applies linear regression to measurements over time to test for a significant degradation trend (ISO Guide 35 classical stability approach).",
        "stab_time_label": "Time Points (e.g. day/month — one per line)",
        "stab_value_label": "Measured Values (paired with time points, same order)",
        "stab_duration_label": "Study Duration (shelf life, same units as time points)",
        "stab_calc_btn": "📊 Test Stability",
        "stab_result_stable": "✅ **Stable** (slope p={p:.4f} ≥ 0.05). No statistically significant trend over time was detected.",
        "stab_result_unstable": "⚠️ **Significant Trend Detected** (slope p={p:.4f} < 0.05). The material may be changing over time — review the shelf-life claim.",
        "stab_ustab_label": "Stability Uncertainty Contribution (u_stab)",
        "stab_slope_label": "Slope (unit/time)",
        "stab_chart_title": "Stability — Value vs Time",
        "unc_description": "Combines characterization, homogeneity, and stability uncertainty components using the GUM approach (sum of squares) to compute the final expanded uncertainty. If you've run the Homogeneity/Stability tools, values are auto-filled.",
        "unc_assigned_value_label": "Assigned Value",
        "unc_u_char_label": "Characterization Uncertainty (u_char, absolute units)",
        "unc_u_char_help": "Standard uncertainty of the measurements used to determine the assigned value (e.g. standard error of the mean).",
        "unc_u_bb_label": "Homogeneity Uncertainty (u_bb)",
        "unc_u_stab_label": "Stability Uncertainty (u_stab)",
        "unc_k_label": "Coverage Factor (k)",
        "unc_use_cached": "Auto-fill from Homogeneity/Stability tabs",
        "unc_calc_btn": "📊 Calculate Uncertainty Budget",
        "unc_result_title": "Assigned Value = {value:.5f} ± {U:.5f} (k={k}, {urel:.2f}% relative)",
        "unc_no_cache": "ℹ️ No previously calculated Homogeneity/Stability result found — manual input will be used.",
        "coa_description": "Generates a formal-style Certificate of Analysis (CoA) PDF containing the assigned value, expanded uncertainty, and traceability information.",
        "coa_material_name": "Material Name",
        "coa_lot_number": "Batch/Lot Number",
        "coa_producer": "Producer/Laboratory",
        "coa_assigned_value_label": "Assigned Value",
        "coa_unit_label": "Unit",
        "coa_expanded_unc_label": "Expanded Uncertainty (U)",
        "coa_k_label": "Coverage Factor (k)",
        "coa_traceability_label": "Traceability Statement",
        "coa_traceability_default": "Traceability of this value is based on absolute quantification via digital PCR using Poisson statistics, traceable to the SI unit (mole).",
        "coa_validity_label": "Expiry Date / Shelf Life",
        "coa_generate_btn": "📥 Generate Certificate (PDF)",
        "coa_download_btn": "⬇️ Certificate of Analysis (CoA)",
        "tab_batch": "Batch Screening",
        "batch_title": "🔬 Batch Sample Screening (CNV Screening)",
        "batch_description": "Designed for quickly screening many samples (e.g. a cohort) against a single reference/expected ratio. Instead of building replicate groups per sample, upload a CSV and a Poisson-based confidence interval is computed for each sample individually.",
        "batch_uploader": "Upload batch screening CSV file",
        "batch_assay_title": "**Assay Assignment**",
        "batch_target_label": "Target gene assay",
        "batch_ref_label": "Reference gene assay(s) (partitions pooled if multiple selected)",
        "batch_settings_title": "**Screening Settings**",
        "batch_expected_ratio": "Expected Ratio (no-change reference)",
        "batch_ploidy": "Reference Ploidy",
        "batch_run_btn": "▶ Run Screening",
        "batch_results_title": "📋 Screening Results",
        "batch_col_sample": "Sample",
        "batch_col_lambda_t": "λ Target",
        "batch_col_lambda_r": "λ Reference",
        "batch_col_ratio": "Ratio",
        "batch_col_ci": "95% CI (Ratio)",
        "batch_col_cn": "Copy Number (CN)",
        "batch_col_conc": "Concentration (copies/µL)",
        "batch_col_class": "Classification",
        "batch_class_gain": "📈 Gain",
        "batch_class_loss": "📉 Loss",
        "batch_class_normal": "➖ Normal",
        "batch_chart_title": "Copy Number Screening Across Samples",
        "batch_download_btn": "📥 Download Screening Results (CSV)",
        "batch_no_data": "No results yet — upload a file and run the screening first.",
        "batch_n_flagged": "{n} / {total} samples are statistically significantly different from the expected ratio (95% CI does not include the no-change value).",
        "tab_vaf": "VAF Calculator",
        "vaf_title": "🧬 VAF / Mutation Fraction Calculator",
        "vaf_description": "Designed for liquid biopsy / ctDNA monitoring. Upload a CSV containing Mutant and Wild-type assays; the Fractional Abundance (FA%, variant allele fraction) with 95% CI is calculated for each sample.",
        "vaf_uploader": "Upload VAF CSV file",
        "vaf_assay_title": "**Assay Assignment**",
        "vaf_mutant_label": "Mutant assay",
        "vaf_wt_label": "Wild-type assay",
        "vaf_run_btn": "▶ Calculate VAF",
        "vaf_results_title": "📋 VAF Results",
        "vaf_col_sample": "Sample",
        "vaf_col_lambda_mut": "λ Mutant",
        "vaf_col_lambda_wt": "λ Wild-type",
        "vaf_col_fa": "FA (%)",
        "vaf_col_ci": "95% CI (FA%)",
        "vaf_col_conc_mut": "Mutant Concentration (copies/µL)",
        "vaf_col_detected": "Detection Status",
        "vaf_detected_yes": "✅ Detected",
        "vaf_detected_no": "❌ Not Detected",
        "vaf_chart_title": "Fractional Abundance (VAF%) Across Samples",
        "vaf_download_btn": "📥 Download VAF Results (CSV)",
        "vaf_no_data": "No results yet — upload a file and run the calculation first.",
        "vaf_method_note": "ℹ️ FA = λ(mutant) / (λ(mutant) + λ(wild-type)). 95% CI calculated via the delta method (Hindson et al. 2013, Anal Chem). This assumes mutant and wild-type assays were measured on the same partition set (multiplex) or under equivalent conditions.",
        "vaf_n_detected": "{n} / {total} samples had the mutant allele detected (≥1 positive partition).",
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

def mean_ci(values):
    """
    Mean, 95% CI (t-distribution, based on replicate-to-replicate variability),
    and CV% for a set of replicate values (e.g. concentration or ratio).
    Returns (mean, ci_low, ci_high, cv_pct, n).
    """
    values = np.array(values, dtype=float)
    values = values[~np.isnan(values)]
    n = len(values)
    if n == 0:
        return np.nan, np.nan, np.nan, np.nan, 0
    m = float(np.mean(values))
    if n < 2:
        return m, np.nan, np.nan, np.nan, n
    sd = float(np.std(values, ddof=1))
    sem = sd / np.sqrt(n)
    t_crit = stats.t.ppf(0.975, df=n - 1)
    ci_low, ci_high = m - t_crit * sem, m + t_crit * sem
    cv_pct = (sd / m * 100) if m != 0 else np.nan
    return m, ci_low, ci_high, cv_pct, n

def lambda_to_conc(lam_array, partition_vol_nl_local):
    """Convert lambda (copies/partition) to concentration in copies/µL."""
    lam_array = np.array(lam_array, dtype=float)
    return lam_array / partition_vol_nl_local * 1000.0

def poisson_ratio_ci(lam_t, se_t, lam_r, se_r):
    """
    95% CI for the ratio R = lambda_t / lambda_r of two independent
    Poisson-derived rates, using the log-scale delta method:
        Var(ln R) ~= (SE_t/lambda_t)^2 + (SE_r/lambda_r)^2
    This avoids negative CI bounds (unlike a linear delta method) and is the
    standard approach for ratios of counting-statistics-derived quantities.
    Returns (ratio, ci_low, ci_high).
    """
    if lam_t is None or lam_r is None or lam_r == 0 or np.isnan(lam_t) or np.isnan(lam_r):
        return np.nan, np.nan, np.nan
    ratio = lam_t / lam_r
    if se_t is None or se_r is None or np.isnan(se_t) or np.isnan(se_r) or lam_t == 0:
        return ratio, np.nan, np.nan
    var_ln_r = (se_t / lam_t) ** 2 + (se_r / lam_r) ** 2
    se_ln_r = np.sqrt(var_ln_r)
    ci_low = ratio * np.exp(-1.96 * se_ln_r)
    ci_high = ratio * np.exp(1.96 * se_ln_r)
    return ratio, ci_low, ci_high

def pool_and_compute_batch(std_df, target_assay, ref_assays, partition_vol_nl_local):
    """
    For each unique Sample in std_df, pools (sums) positive/total partition
    counts across replicate rows for the target assay and for the combined
    set of reference assays, then computes lambda/SE for each and the
    Ratio with its 95% CI via poisson_ratio_ci.
    Returns a list of per-sample result dicts.
    """
    results = []
    samples = sorted(std_df["Sample"].unique())
    for sample in samples:
        tgt_rows = std_df[(std_df["Sample"] == sample) & (std_df["Target"] == target_assay)]
        ref_rows = std_df[(std_df["Sample"] == sample) & (std_df["Target"].isin(ref_assays))]
        if tgt_rows.empty or ref_rows.empty:
            continue

        pos_t, tot_t = float(tgt_rows["Positives"].sum()), float(tgt_rows["Total"].sum())
        pos_r, tot_r = float(ref_rows["Positives"].sum()), float(ref_rows["Total"].sum())

        lam_t, _, _, status_t = poisson_lambda(pos_t, tot_t)
        lam_r, _, _, status_r = poisson_lambda(pos_r, tot_r)
        if status_t != "ok" or status_r != "ok":
            continue

        p_t, p_r = pos_t / tot_t, pos_r / tot_r
        se_t = np.sqrt(p_t / (tot_t * (1 - p_t))) if (1 - p_t) > 0 else np.nan
        se_r = np.sqrt(p_r / (tot_r * (1 - p_r))) if (1 - p_r) > 0 else np.nan

        ratio, ci_low, ci_high = poisson_ratio_ci(lam_t, se_t, lam_r, se_r)
        conc = lam_t / partition_vol_nl_local * 1000.0

        results.append({
            "Sample": sample, "lambda_t": lam_t, "lambda_r": lam_r,
            "ratio": ratio, "ci_low": ci_low, "ci_high": ci_high, "conc": conc,
            "n_partitions_t": int(tot_t), "n_partitions_r": int(tot_r),
        })
    return results

def compute_vaf(lam_mut, se_mut, lam_wt, se_wt):
    """
    Fractional Abundance (FA, i.e. variant allele fraction) for a mutant vs
    wild-type ddPCR assay pair: FA = lambda_mut / (lambda_mut + lambda_wt).
    95% CI via the linear delta method on the ratio-to-sum form, which is
    the standard approach used in ddPCR liquid-biopsy literature
    (e.g. Hindson et al. 2013 Anal Chem):
        Var(FA) ~= [y/(x+y)^2]^2 * Var(x) + [x/(x+y)^2]^2 * Var(y)
    Returns (FA, ci_low, ci_high) as fractions in [0, 1].
    """
    x, y = lam_mut, lam_wt
    if x is None or y is None or np.isnan(x) or np.isnan(y) or (x + y) <= 0:
        return np.nan, np.nan, np.nan
    fa = x / (x + y)
    if se_mut is None or se_wt is None or np.isnan(se_mut) or np.isnan(se_wt):
        return fa, np.nan, np.nan
    denom_sq = (x + y) ** 4
    var_fa = ((y ** 2) * (se_mut ** 2) + (x ** 2) * (se_wt ** 2)) / denom_sq
    se_fa = np.sqrt(var_fa)
    ci_low = max(fa - 1.96 * se_fa, 0.0)
    ci_high = min(fa + 1.96 * se_fa, 1.0)
    return fa, ci_low, ci_high

def pool_and_compute_vaf(std_df, mutant_assay, wt_assay, partition_vol_nl_local):
    """
    For each unique Sample in std_df, pools positive/total partition counts
    for the mutant assay and the wild-type assay, computes lambda/SE for
    each, and the Fractional Abundance (FA%) with 95% CI via compute_vaf.
    Returns a list of per-sample result dicts.
    """
    results = []
    samples = sorted(std_df["Sample"].unique())
    for sample in samples:
        mut_rows = std_df[(std_df["Sample"] == sample) & (std_df["Target"] == mutant_assay)]
        wt_rows = std_df[(std_df["Sample"] == sample) & (std_df["Target"] == wt_assay)]
        if mut_rows.empty or wt_rows.empty:
            continue

        pos_m, tot_m = float(mut_rows["Positives"].sum()), float(mut_rows["Total"].sum())
        pos_w, tot_w = float(wt_rows["Positives"].sum()), float(wt_rows["Total"].sum())

        lam_m, _, _, status_m = poisson_lambda(pos_m, tot_m)
        lam_w, _, _, status_w = poisson_lambda(pos_w, tot_w)
        if status_m not in ("ok", "saturated") or status_w not in ("ok", "saturated"):
            continue
        if status_m == "saturated":
            lam_m = -np.log(1.0 / tot_m)  # fallback, extreme edge case
        if status_w == "saturated":
            lam_w = -np.log(1.0 / tot_w)
        lam_m = 0.0 if pos_m == 0 else lam_m
        lam_w = 0.0 if pos_w == 0 else lam_w

        p_m = pos_m / tot_m if tot_m > 0 else np.nan
        p_w = pos_w / tot_w if tot_w > 0 else np.nan
        se_m = np.sqrt(p_m / (tot_m * (1 - p_m))) if (1 - p_m) > 0 else np.nan
        se_w = np.sqrt(p_w / (tot_w * (1 - p_w))) if (1 - p_w) > 0 else np.nan

        fa, ci_low, ci_high = compute_vaf(lam_m, se_m, lam_w, se_w)
        conc_mut = lam_m / partition_vol_nl_local * 1000.0
        detected = pos_m >= 1

        results.append({
            "Sample": sample, "lambda_mut": lam_m, "lambda_wt": lam_w,
            "fa": fa, "ci_low": ci_low, "ci_high": ci_high,
            "conc_mut": conc_mut, "detected": detected,
            "n_partitions_mut": int(tot_m), "n_partitions_wt": int(tot_w),
            "pos_mut": int(pos_m),
        })
    return results

_PROJECT_SCALAR_KEYS = [
    "gene_count", "patient_count", "num_ref_genes", "ploidy", "partition_vol",
    "qc_min", "outlier_enabled", "outlier_method", "grubbs_alpha", "iqr_mult",
]
_PROJECT_TEXT_PREFIXES = (
    "ctrl_tgt_pos_", "ctrl_tgt_tot_", "ctrl_ref_pos_", "ctrl_ref_tot_",
    "smp_tgt_pos_", "smp_tgt_tot_", "smp_ref_pos_", "smp_ref_tot_",
    "ntc_pos_", "ntc_tot_",
)

def export_project_state():
    """Serializes all data-entry-relevant session_state keys into a JSON-safe dict."""
    project = {"_absolutegene_project_version": 1}
    for k in _PROJECT_SCALAR_KEYS:
        if k in st.session_state:
            v = st.session_state[k]
            project[k] = v if isinstance(v, (int, float, str, bool)) else str(v)
    for k in list(st.session_state.keys()):
        if isinstance(k, str) and k.startswith(_PROJECT_TEXT_PREFIXES):
            v = st.session_state[k]
            if isinstance(v, str):
                project[k] = v
    return project

def import_project_state(project_dict):
    """Restores session_state from a previously exported project dict. Returns count of keys restored."""
    if not isinstance(project_dict, dict):
        return 0
    count = 0
    for k, v in project_dict.items():
        if k == "_absolutegene_project_version":
            continue
        if k in _PROJECT_SCALAR_KEYS or (isinstance(k, str) and k.startswith(_PROJECT_TEXT_PREFIXES)):
            st.session_state[k] = v
            count += 1
    return count

def recommend_dilution(current_lambda, target_lambda=1.6):
    """
    Given the current observed lambda (copies/partition), recommends a
    dilution adjustment factor to bring the assay into the optimal dynamic
    range for Poisson-based quantification.

    Ideal operating range: lambda ~0.05-3.0 (acceptable), ~1.0-2.0 (optimal
    precision/dynamic-range balance; default target 1.6, near the point of
    minimal relative CI width for Poisson-derived concentration estimates).

    Returns dict: status ('too_low'/'optimal'/'too_high'/'saturated_risk'),
    dilution_factor (multiply current dilution by this value to reach
    target_lambda; >1 means dilute further, <1 means concentrate/dilute less),
    message_key (translation key hint for status).
    """
    if current_lambda is None or np.isnan(current_lambda) or current_lambda <= 0:
        return None
    dilution_factor = current_lambda / target_lambda
    if current_lambda < 0.05:
        status = "too_low"
    elif current_lambda > 4.0:
        status = "saturated_risk"
    elif current_lambda > 3.0:
        status = "too_high"
    else:
        status = "optimal"
    return {
        "current_lambda": current_lambda, "target_lambda": target_lambda,
        "dilution_factor": dilution_factor, "status": status,
    }

# ═══════════════════════════════════════════════════════════════════════════════
# CLINICAL TOOLS — Measurement Uncertainty, RCV, Precision Study, Method Comparison
# ═══════════════════════════════════════════════════════════════════════════════
def compute_mu_budget(u_poisson_pct, u_pipetting_pct, u_precision_pct, k=2.0):
    """
    Combined and expanded measurement uncertainty following the GUM
    (Guide to the Expression of Uncertainty in Measurement) approach:
    relative standard uncertainties from independent sources are combined
    in quadrature, then multiplied by a coverage factor k (k=2 for an
    approximate 95% confidence level, assuming a normal distribution).
        u_c(%) = sqrt(u_poisson^2 + u_pipetting^2 + u_precision^2)
        U(%) = k * u_c
    Returns dict with u_c_pct, U_pct, and the individual contributions
    (useful for a budget table showing each source's relative share).
    """
    components = {
        "Poisson counting statistics": u_poisson_pct or 0.0,
        "Pipetting / dilution": u_pipetting_pct or 0.0,
        "Inter-run precision": u_precision_pct or 0.0,
    }
    u_c_pct = np.sqrt(sum(v ** 2 for v in components.values()))
    U_pct = k * u_c_pct
    return {"components": components, "u_c_pct": u_c_pct, "U_pct": U_pct, "k": k}

def compute_rcv(result1, result2, cv_analytical_pct, cv_biological_pct, z=1.96):
    """
    Reference Change Value (RCV) — determines whether the difference between
    two serial results exceeds what would be expected from analytical and
    within-subject biological variation alone (Fraser CG, Harris EK.
    Generation and application of data on biological variation in clinical
    chemistry. Crit Rev Clin Lab Sci 1989).
        RCV(%) = sqrt(2) * z * sqrt(CV_analytical^2 + CV_biological^2)
    A percent change between result1 and result2 exceeding RCV% is
    considered a statistically significant change beyond combined
    analytical + biological noise.
    Returns dict with rcv_pct, percent_change, significant (bool).
    """
    if result1 is None or result2 is None or result1 == 0:
        return None
    rcv_pct = np.sqrt(2) * z * np.sqrt(cv_analytical_pct ** 2 + cv_biological_pct ** 2)
    percent_change = (result2 - result1) / abs(result1) * 100.0
    significant = abs(percent_change) > rcv_pct
    return {"rcv_pct": rcv_pct, "percent_change": percent_change, "significant": significant}

def compute_precision_study(day_groups):
    """
    Nested one-way ANOVA (day/run as the grouping factor) for a basic
    precision (repeatability / intermediate precision) study, following
    the general approach of CLSI EP05-A3. Assumes equal replicates per day
    for simplicity (unbalanced designs are common in practice but require
    more elaborate variance-component estimation).

    day_groups: list of arrays, one array of replicate values per day/run.

    Returns dict with repeatability_cv, between_day_cv, total_cv (all %),
    grand_mean, and the underlying ANOVA mean squares.
    """
    day_groups = [np.array(g, dtype=float) for g in day_groups if len(g) > 0]
    k = len(day_groups)
    if k < 2:
        return None
    n_per_day = [len(g) for g in day_groups]
    if len(set(n_per_day)) != 1:
        n = min(n_per_day)
        day_groups = [g[:n] for g in day_groups]
    else:
        n = n_per_day[0]
    if n < 2:
        return None

    all_values = np.concatenate(day_groups)
    grand_mean = np.mean(all_values)
    day_means = np.array([np.mean(g) for g in day_groups])

    ss_within = sum(np.sum((g - np.mean(g)) ** 2) for g in day_groups)
    df_within = k * (n - 1)
    ms_within = ss_within / df_within if df_within > 0 else np.nan

    ss_between = n * np.sum((day_means - grand_mean) ** 2)
    df_between = k - 1
    ms_between = ss_between / df_between if df_between > 0 else np.nan

    s_between_day_sq = max((ms_between - ms_within) / n, 0.0) if not np.isnan(ms_between) else 0.0
    s_repeatability = np.sqrt(ms_within) if ms_within >= 0 else np.nan
    s_total = np.sqrt(ms_within + s_between_day_sq)

    return {
        "grand_mean": grand_mean, "k_days": k, "n_per_day": n,
        "ms_within": ms_within, "ms_between": ms_between,
        "repeatability_cv": (s_repeatability / grand_mean * 100) if grand_mean != 0 else np.nan,
        "between_day_cv": (np.sqrt(s_between_day_sq) / grand_mean * 100) if grand_mean != 0 else np.nan,
        "total_cv": (s_total / grand_mean * 100) if grand_mean != 0 else np.nan,
    }

def compute_bland_altman(method1_vals, method2_vals):
    """
    Bland-Altman agreement analysis for paired measurements from two
    methods (Bland JM, Altman DG. Statistical methods for assessing
    agreement between two methods of clinical measurement. Lancet 1986).
    Returns dict with mean_diff (bias), sd_diff, loa_low, loa_high,
    means (array of per-pair means), diffs (array of per-pair differences).
    """
    m1 = np.array(method1_vals, dtype=float)
    m2 = np.array(method2_vals, dtype=float)
    n = min(len(m1), len(m2))
    m1, m2 = m1[:n], m2[:n]
    if n < 2:
        return None
    diffs = m2 - m1
    means = (m1 + m2) / 2
    mean_diff = float(np.mean(diffs))
    sd_diff = float(np.std(diffs, ddof=1))
    return {
        "means": means, "diffs": diffs, "mean_diff": mean_diff, "sd_diff": sd_diff,
        "loa_low": mean_diff - 1.96 * sd_diff, "loa_high": mean_diff + 1.96 * sd_diff,
        "n": n,
    }

def compute_deming_regression(method1_vals, method2_vals, variance_ratio=1.0):
    """
    Deming regression: linear regression accounting for measurement error
    in both variables (unlike ordinary least squares, which assumes the
    x-variable is error-free). variance_ratio (lambda) = Var(error_y)/Var(error_x);
    lambda=1 assumes equal error variance in both methods (a common default
    when no independent estimate of the error ratio is available).
    Returns dict with slope, intercept.
    Reference: Linnet K. Estimation of the linear relationship between the
    measurements of two methods with proportional errors. Stat Med 1990.
    """
    x = np.array(method1_vals, dtype=float)
    y = np.array(method2_vals, dtype=float)
    n = min(len(x), len(y))
    x, y = x[:n], y[:n]
    if n < 3:
        return None
    mx, my = np.mean(x), np.mean(y)
    sxx = np.sum((x - mx) ** 2) / (n - 1)
    syy = np.sum((y - my) ** 2) / (n - 1)
    sxy = np.sum((x - mx) * (y - my)) / (n - 1)
    lam = variance_ratio
    slope = (syy - lam * sxx + np.sqrt((syy - lam * sxx) ** 2 + 4 * lam * sxy ** 2)) / (2 * sxy) if sxy != 0 else np.nan
    intercept = my - slope * mx
    return {"slope": slope, "intercept": intercept, "n": n}

# ═══════════════════════════════════════════════════════════════════════════════
# CRM PRODUCTION — Homogeneity, Stability, Assigned Value/Uncertainty Budget (per
# ISO Guide 35 / ISO 17034 general principles), and Certificate generation.
# ═══════════════════════════════════════════════════════════════════════════════
def compute_homogeneity(unit_groups):
    """
    One-way ANOVA-based between-unit homogeneity test, following the general
    approach described in ISO Guide 35 and Linsinger TP, Pauwels J, van der
    Veen AMH, Schimmel H, Lamberty A. "Homogeneity and stability of reference
    materials." Accred Qual Assur 2001;6:20-25.

    unit_groups: list of arrays, one array of replicate measurements per
    CRM unit/vial (equal replicate count per unit assumed for simplicity).

    Returns dict with ms_within, ms_between, F, F_crit, p_value,
    is_homogeneous (bool, F <= F_crit at alpha=0.05), s_bb (between-unit
    standard deviation estimate), u_bb (standard uncertainty contribution
    from potential inhomogeneity), grand_mean.
    """
    unit_groups = [np.array(g, dtype=float) for g in unit_groups if len(g) > 0]
    p = len(unit_groups)
    if p < 2:
        return None
    n_per_unit = [len(g) for g in unit_groups]
    if len(set(n_per_unit)) != 1:
        n = min(n_per_unit)
        unit_groups = [g[:n] for g in unit_groups]
    else:
        n = n_per_unit[0]
    if n < 2:
        return None

    all_values = np.concatenate(unit_groups)
    grand_mean = float(np.mean(all_values))
    unit_means = np.array([np.mean(g) for g in unit_groups])

    ss_within = sum(np.sum((g - np.mean(g)) ** 2) for g in unit_groups)
    df_within = p * (n - 1)
    ms_within = ss_within / df_within if df_within > 0 else np.nan

    ss_between = n * np.sum((unit_means - grand_mean) ** 2)
    df_between = p - 1
    ms_between = ss_between / df_between if df_between > 0 else np.nan

    F = ms_between / ms_within if (ms_within and ms_within > 0) else np.nan
    F_crit = stats.f.ppf(0.95, df_between, df_within) if df_within > 0 else np.nan
    p_value = stats.f.sf(F, df_between, df_within) if not np.isnan(F) else np.nan
    is_homogeneous = (F <= F_crit) if not np.isnan(F) else True

    if ms_between > ms_within:
        s_bb_sq = (ms_between - ms_within) / n
        u_bb = np.sqrt(s_bb_sq)
    else:
        # Conservative minimum-uncertainty estimate when between-unit
        # variance is not resolvable from within-unit noise (Linsinger et al. 2001)
        u_bb = np.sqrt(ms_within / n) * (2.0 / df_between) ** 0.25 if df_between > 0 else np.nan
        s_bb_sq = 0.0

    return {
        "grand_mean": grand_mean, "p_units": p, "n_per_unit": n,
        "ms_within": ms_within, "ms_between": ms_between, "F": F, "F_crit": F_crit,
        "p_value": p_value, "is_homogeneous": is_homogeneous,
        "s_bb": np.sqrt(s_bb_sq), "u_bb": u_bb,
        "df_within": df_within, "df_between": df_between,
    }

def compute_stability(time_points, values, study_duration=None):
    """
    Linear-regression-based stability assessment: tests whether the analyte
    concentration/ratio shows a statistically significant trend over time
    (ISO Guide 35 classical stability study approach). If a significant
    trend is found, this may indicate degradation and should inform the
    material's shelf-life claim.

    time_points, values: paired arrays (e.g. storage duration in days/months
    and the corresponding measured value).
    study_duration: the shelf-life/study period of interest, used to compute
    the uncertainty contribution from potential instability
    (u_stab = SE(slope) * study_duration). If None, uses max(time_points).

    Returns dict with slope, intercept, se_slope, t_stat, p_value,
    is_stable (bool, p >= 0.05), u_stab.
    """
    x = np.array(time_points, dtype=float)
    y = np.array(values, dtype=float)
    n = min(len(x), len(y))
    x, y = x[:n], y[:n]
    if n < 3:
        return None

    slope, intercept, r_value, p_value, se_slope = stats.linregress(x, y)
    t_stat = slope / se_slope if se_slope > 0 else np.nan
    is_stable = p_value >= 0.05

    duration = study_duration if study_duration is not None else float(np.max(x))
    u_stab = abs(se_slope) * duration

    return {
        "slope": slope, "intercept": intercept, "se_slope": se_slope,
        "t_stat": t_stat, "p_value": p_value, "r_value": r_value,
        "is_stable": is_stable, "u_stab": u_stab, "study_duration": duration, "n": n,
    }

def compute_assigned_value_uncertainty(assigned_value, u_char, u_bb, u_stab,
                                        extra_components=None, k=2.0):
    """
    Combines characterization uncertainty (u_char), between-unit
    inhomogeneity uncertainty (u_bb), stability uncertainty (u_stab), and
    any additional user-specified components into a combined standard
    uncertainty and expanded uncertainty, following the GUM approach
    (components combined in quadrature, assuming independence):
        u_c = sqrt(u_char^2 + u_bb^2 + u_stab^2 + sum(extra_i^2))
        U = k * u_c
    All uncertainty inputs are absolute (same units as assigned_value), not
    relative percentages, consistent with ISO Guide 35 uncertainty budgets
    for certified reference materials.
    Returns dict with components, u_c, U, U_rel_pct.
    """
    components = {"Characterization": u_char or 0.0, "Homogeneity (u_bb)": u_bb or 0.0,
                  "Stability (u_stab)": u_stab or 0.0}
    if extra_components:
        components.update(extra_components)
    u_c = np.sqrt(sum(v ** 2 for v in components.values()))
    U = k * u_c
    U_rel_pct = (U / assigned_value * 100) if assigned_value else np.nan
    return {"components": components, "u_c": u_c, "U": U, "U_rel_pct": U_rel_pct, "k": k,
            "assigned_value": assigned_value}

def compute_lod_loq(ntc_pos, ntc_tot, partition_vol_nl_local):
    """
    Estimate LOD (limit of detection) and LOQ (limit of quantification) from
    pooled No-Template Control (NTC) replicate partition counts.

    If zero positives are observed across all pooled NTC partitions, LOD is
    derived from the "rule of three" upper 95% confidence bound on the true
    positive rate: p_upper ≈ 3/n_total. If NTC shows background positives
    (possible contamination), LOD is set to the NTC's own λ plus a one-sided
    95% margin (background + 1.645×SEM).

    LOQ is reported as 3×LOD, a commonly used simplified heuristic; true LOQ
    should ideally be validated empirically via replicate reproducibility (CV%).

    Returns dict with keys: lod_lambda, loq_lambda, lod_conc, loq_conc,
    ntc_lambda, ntc_pos_total, ntc_tot_total, contamination (bool).
    """
    ntc_pos = np.array(ntc_pos, dtype=float)
    ntc_tot = np.array(ntc_tot, dtype=float)
    n = min(len(ntc_pos), len(ntc_tot))
    if n == 0:
        return None
    ntc_pos, ntc_tot = ntc_pos[:n], ntc_tot[:n]

    pos_total = float(np.sum(ntc_pos))
    tot_total = float(np.sum(ntc_tot))
    if tot_total <= 0:
        return None

    if pos_total == 0:
        # Rule-of-three upper bound on the true positive rate
        p_upper = 3.0 / tot_total
        lod_lambda = -np.log(1 - p_upper) if p_upper < 1 else np.nan
        ntc_lambda = 0.0
        contamination = False
    else:
        p = pos_total / tot_total
        ntc_lambda, _, _, _ = poisson_lambda(pos_total, tot_total)
        se = np.sqrt(p / (tot_total * (1 - p))) if (1 - p) > 0 else np.nan
        lod_lambda = ntc_lambda + 1.645 * se if not np.isnan(se) else ntc_lambda
        contamination = True

    loq_lambda = lod_lambda * 3.0
    return {
        "lod_lambda": lod_lambda, "loq_lambda": loq_lambda,
        "lod_conc": lod_lambda / partition_vol_nl_local * 1000.0,
        "loq_conc": loq_lambda / partition_vol_nl_local * 1000.0,
        "ntc_lambda": ntc_lambda, "ntc_pos_total": pos_total, "ntc_tot_total": tot_total,
        "contamination": contamination,
    }

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
# INSTRUMENT CSV IMPORT (QuantaSoft / QX Manager / QIAcuity / generic)
# ═══════════════════════════════════════════════════════════════════════════════
_COLUMN_ALIASES = {
    "sample": ["sample", "samplename", "sample name", "well name", "wellname",
               "sample description 1", "sample_id", "sampleid"],
    "target": ["target", "target name", "targettype", "target type", "assay",
               "biomarker", "assay name", "gene", "channel"],
    "positives": ["positives", "positive droplets", "positives (rain excluded)",
                  "ch1 positives", "ch1+", "positivecount", "positive count",
                  "positive partitions", "positive reactions"],
    "total": ["accepteddroplets", "accepted droplets", "total droplets",
              "total accepted droplets", "valid partitions", "partitions (valid)",
              "total", "total partitions", "total reactions", "total conf. droplets",
              "droplet count", "well droplet count"],
}

def _norm_col(col):
    return str(col).strip().lower()

def _auto_detect_column(df_columns, field):
    normed = {_norm_col(c): c for c in df_columns}
    for alias in _COLUMN_ALIASES[field]:
        if alias in normed:
            return normed[alias]
    # partial/substring match fallback
    for norm_c, orig_c in normed.items():
        for alias in _COLUMN_ALIASES[field]:
            if alias in norm_c or norm_c in alias:
                return orig_c
    return None

def parse_instrument_csv(file_bytes):
    """
    Parses a generic dPCR/ddPCR instrument export (QuantaSoft, QX Manager,
    QIAcuity, or similar CSV/TSV exports). Auto-detects Sample, Target,
    Positives, and Total (accepted partitions) columns using common column
    name aliases across instruments. Returns (df, detected_cols, error) where
    df has standardized columns: Sample, Target, Positives, Total.
    detected_cols is a dict {field: column_name_or_None} for any fields that
    need manual mapping.
    """
    import io as _io
    try:
        content = file_bytes.decode("utf-8-sig", errors="replace")
    except Exception as e:
        return None, None, f"File decoding error: {e}"

    sep = "\t" if content.count("\t") > content.count(",") else ","
    try:
        df = pd.read_csv(_io.StringIO(content), sep=sep)
    except Exception as e:
        return None, None, f"CSV parse error: {e}"

    if df.empty or len(df.columns) < 2:
        return None, None, "No usable columns found in file."

    detected = {field: _auto_detect_column(df.columns, field) for field in _COLUMN_ALIASES}
    return df, detected, None

def build_standard_import_df(raw_df, col_map):
    """
    col_map: {"sample": colname, "target": colname, "positives": colname, "total": colname}
    Returns a cleaned DataFrame with columns Sample, Target, Positives, Total (numeric).
    """
    out = pd.DataFrame()
    out["Sample"] = raw_df[col_map["sample"]].astype(str).str.strip()
    out["Target"] = raw_df[col_map["target"]].astype(str).str.strip()
    out["Positives"] = pd.to_numeric(raw_df[col_map["positives"]], errors="coerce")
    out["Total"] = pd.to_numeric(raw_df[col_map["total"]], errors="coerce")
    out = out.dropna(subset=["Positives", "Total"])
    out = out[out["Total"] > 0]
    return out.reset_index(drop=True)

def apply_csv_import_to_session(std_df, target_assays, ref_assays, ctrl_label, patient_labels):
    """
    Given the standardized import DataFrame (Sample, Target, Positives, Total)
    and the user's assay/group assignments, fills st.session_state text_area
    values for each target gene / reference gene / control / patient group
    combination. Returns the number of individual partition-count values filled.
    """
    if std_df is None or std_df.empty:
        return 0

    ctrl_keywords = [s.strip() for s in ctrl_label.split(",") if s.strip()]

    def is_ctrl(name):
        return any(kw.lower() in name.lower() for kw in ctrl_keywords)

    count = 0
    for gi, target_name in enumerate(target_assays):
        tgt_df = std_df[std_df["Target"] == target_name]

        ctrl_rows = tgt_df[tgt_df["Sample"].apply(is_ctrl)]
        if len(ctrl_rows) > 0:
            st.session_state[f"ctrl_tgt_pos_{gi}"] = "\n".join(str(int(v)) for v in ctrl_rows["Positives"])
            st.session_state[f"ctrl_tgt_tot_{gi}"] = "\n".join(str(int(v)) for v in ctrl_rows["Total"])
            count += len(ctrl_rows) * 2

        for ri, ref_name in enumerate(ref_assays):
            ref_df = std_df[std_df["Target"] == ref_name]
            ref_ctrl_rows = ref_df[ref_df["Sample"].apply(is_ctrl)]
            if len(ref_ctrl_rows) > 0:
                st.session_state[f"ctrl_ref_pos_{gi}_{ri}"] = "\n".join(str(int(v)) for v in ref_ctrl_rows["Positives"])
                st.session_state[f"ctrl_ref_tot_{gi}_{ri}"] = "\n".join(str(int(v)) for v in ref_ctrl_rows["Total"])
                count += len(ref_ctrl_rows) * 2

        for pj, pat_label in enumerate(patient_labels):
            pat_keywords = [s.strip() for s in pat_label.split(",") if s.strip()]

            def is_pat(name, kws=pat_keywords):
                return any(kw.lower() in name.lower() for kw in kws)

            pat_rows = tgt_df[tgt_df["Sample"].apply(is_pat)]
            if len(pat_rows) > 0:
                st.session_state[f"smp_tgt_pos_{gi}_{pj}"] = "\n".join(str(int(v)) for v in pat_rows["Positives"])
                st.session_state[f"smp_tgt_tot_{gi}_{pj}"] = "\n".join(str(int(v)) for v in pat_rows["Total"])
                count += len(pat_rows) * 2

            for ri, ref_name in enumerate(ref_assays):
                ref_df = std_df[std_df["Target"] == ref_name]
                ref_pat_rows = ref_df[ref_df["Sample"].apply(is_pat)]
                if len(ref_pat_rows) > 0:
                    st.session_state[f"smp_ref_pos_{gi}_{pj}_{ri}"] = "\n".join(str(int(v)) for v in ref_pat_rows["Positives"])
                    st.session_state[f"smp_ref_tot_{gi}_{pj}_{ri}"] = "\n".join(str(int(v)) for v in ref_pat_rows["Total"])
                    count += len(ref_pat_rows) * 2

    st.session_state["gene_count"] = max(1, len(target_assays))
    st.session_state["num_ref_genes"] = max(1, len(ref_assays))
    st.session_state["patient_count"] = max(1, len(patient_labels))
    return count

# ═══════════════════════════════════════════════════════════════════════════════
# EXAMPLE CSV FILES (for the three upload-based workflows)
# ═══════════════════════════════════════════════════════════════════════════════
EXAMPLE_CSV_MAIN_IMPORT = """Sample,Target,Positives,AcceptedDroplets
Control_1,MYCN,1890,20000
Control_1,RPP30,1895,20000
Control_2,MYCN,1920,20100
Control_2,RPP30,1915,20100
Control_3,MYCN,1875,19850
Control_3,RPP30,1870,19850
Control_4,MYCN,1905,20050
Control_4,RPP30,1900,20050
Control_5,MYCN,1898,19980
Control_5,RPP30,1892,19980
Patient_1,MYCN,2780,20000
Patient_1,RPP30,1895,20000
Patient_2,MYCN,2820,20100
Patient_2,RPP30,1915,20100
Patient_3,MYCN,2755,19850
Patient_3,RPP30,1870,19850
Patient_4,MYCN,2800,20050
Patient_4,RPP30,1900,20050
Patient_5,MYCN,2790,19980
Patient_5,RPP30,1892,19980
"""

EXAMPLE_CSV_BATCH_SCREENING = """Sample,Target,Positives,AcceptedDroplets
Sample_01,MYCN,1902,20000
Sample_01,RPP30,1898,20000
Sample_02,MYCN,1875,19950
Sample_02,RPP30,1910,19950
Sample_03,MYCN,1940,20100
Sample_03,RPP30,1885,20100
Sample_04,MYCN,1888,19870
Sample_04,RPP30,1901,19870
Sample_05,MYCN,1915,20050
Sample_05,RPP30,1893,20050
Sample_06,MYCN,1867,19920
Sample_06,RPP30,1878,19920
Sample_07,MYCN,1933,20080
Sample_07,RPP30,1922,20080
Sample_08,MYCN,4980,20000
Sample_08,RPP30,1899,20000
Sample_09,MYCN,5120,20100
Sample_09,RPP30,1907,20100
Sample_10,MYCN,4870,19950
Sample_10,RPP30,1884,19950
"""

EXAMPLE_CSV_VAF = """Sample,Target,Positives,AcceptedDroplets
Baseline,Mutant_KRASG12D,450,18000
Baseline,WT_KRAS,1200,18000
Post-Treatment,Mutant_KRASG12D,2,19500
Post-Treatment,WT_KRAS,1350,19500
Follow-up-3mo,Mutant_KRASG12D,0,19800
Follow-up-3mo,WT_KRAS,1400,19800
Relapse,Mutant_KRASG12D,180,18500
Relapse,WT_KRAS,1250,18500
"""

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
# EXAMPLE / VALIDATION SCENARIOS (dPCR partition counts — Poisson-derived)
# ═══════════════════════════════════════════════════════════════════════════════
SCENARIOS = {
    "S1 — Basic CNV gain (1 gene, n=5)": {
        "gene_count": 1, "patient_count": 1, "num_ref_genes": 1,
        "ploidy": 2, "partition_vol": 0.85, "qc_min": 10000,
        "outlier_method": "Grubbs", "outlier_enabled": True,
        "description_tr": "1 hedef gen, 1 hasta grubu, n=5. Hasta grubunda hedef lokusta ~1.5x kopya kazanımı (CN ~3).",
        "description_en": "1 target gene, 1 patient group, n=5. Patient group shows ~1.5x copy gain at the target locus (CN ~3).",
        "ctrl_tgt_pos_0": "1890\n1920\n1875\n1905\n1898",
        "ctrl_tgt_tot_0": "20000\n20100\n19850\n20050\n19980",
        "ctrl_ref_pos_0_0": "1895\n1915\n1870\n1900\n1892",
        "ctrl_ref_tot_0_0": "20000\n20100\n19850\n20050\n19980",
        "smp_tgt_pos_0_0": "2780\n2820\n2755\n2800\n2790",
        "smp_tgt_tot_0_0": "20000\n20100\n19850\n20050\n19980",
        "smp_ref_pos_0_0_0": "1895\n1915\n1870\n1900\n1892",
        "smp_ref_tot_0_0_0": "20000\n20100\n19850\n20050\n19980",
    },
    "S2 — Multi-gene + dual reference (2 genes, 2 groups, n=5)": {
        "gene_count": 2, "patient_count": 2, "num_ref_genes": 2,
        "ploidy": 2, "partition_vol": 0.85, "qc_min": 10000,
        "outlier_method": "IQR", "outlier_enabled": True,
        "description_tr": "2 hedef gen, 2 hasta grubu, n=5. İkili referans normalizasyonu (geNorm). Gen 1: değişim yok. Gen 2: Grup 1 kayıp, Grup 2 kazanım.",
        "description_en": "2 target genes, 2 patient groups, n=5. Dual-reference normalization (geNorm-style). Gene 1: no change. Gene 2: Group 1 = loss, Group 2 = gain.",
        # Gene 1 — control
        "ctrl_tgt_pos_0": "2000\n2020\n1980\n2010\n1995",
        "ctrl_tgt_tot_0": "20000\n20100\n19850\n20050\n19980",
        "ctrl_ref_pos_0_0": "1900\n1920\n1880\n1910\n1898",
        "ctrl_ref_tot_0_0": "20000\n20100\n19850\n20050\n19980",
        "ctrl_ref_pos_0_1": "2100\n2120\n2080\n2110\n2098",
        "ctrl_ref_tot_0_1": "20000\n20100\n19850\n20050\n19980",
        # Gene 1 — Group 1 (no change)
        "smp_tgt_pos_0_0": "2005\n2025\n1985\n2015\n2000",
        "smp_tgt_tot_0_0": "20000\n20100\n19850\n20050\n19980",
        "smp_ref_pos_0_0_0": "1900\n1920\n1880\n1910\n1898",
        "smp_ref_tot_0_0_0": "20000\n20100\n19850\n20050\n19980",
        "smp_ref_pos_0_0_1": "2100\n2120\n2080\n2110\n2098",
        "smp_ref_tot_0_0_1": "20000\n20100\n19850\n20050\n19980",
        # Gene 1 — Group 2 (no change)
        "smp_tgt_pos_0_1": "1995\n2015\n1975\n2005\n1990",
        "smp_tgt_tot_0_1": "20000\n20100\n19850\n20050\n19980",
        "smp_ref_pos_0_1_0": "1900\n1920\n1880\n1910\n1898",
        "smp_ref_tot_0_1_0": "20000\n20100\n19850\n20050\n19980",
        "smp_ref_pos_0_1_1": "2100\n2120\n2080\n2110\n2098",
        "smp_ref_tot_0_1_1": "20000\n20100\n19850\n20050\n19980",
        # Gene 2 — control
        "ctrl_tgt_pos_1": "2000\n2020\n1980\n2010\n1995",
        "ctrl_tgt_tot_1": "20000\n20100\n19850\n20050\n19980",
        "ctrl_ref_pos_1_0": "1900\n1920\n1880\n1910\n1898",
        "ctrl_ref_tot_1_0": "20000\n20100\n19850\n20050\n19980",
        "ctrl_ref_pos_1_1": "2100\n2120\n2080\n2110\n2098",
        "ctrl_ref_tot_1_1": "20000\n20100\n19850\n20050\n19980",
        # Gene 2 — Group 1 (loss, ratio ~0.5)
        "smp_tgt_pos_1_0": "1030\n1045\n1015\n1038\n1022",
        "smp_tgt_tot_1_0": "20000\n20100\n19850\n20050\n19980",
        "smp_ref_pos_1_0_0": "1900\n1920\n1880\n1910\n1898",
        "smp_ref_tot_1_0_0": "20000\n20100\n19850\n20050\n19980",
        "smp_ref_pos_1_0_1": "2100\n2120\n2080\n2110\n2098",
        "smp_ref_tot_1_0_1": "20000\n20100\n19850\n20050\n19980",
        # Gene 2 — Group 2 (gain, ratio ~2.0)
        "smp_tgt_pos_1_1": "3790\n3830\n3760\n3810\n3795",
        "smp_tgt_tot_1_1": "20000\n20100\n19850\n20050\n19980",
        "smp_ref_pos_1_1_0": "1900\n1920\n1880\n1910\n1898",
        "smp_ref_tot_1_1_0": "20000\n20100\n19850\n20050\n19980",
        "smp_ref_pos_1_1_1": "2100\n2120\n2080\n2110\n2098",
        "smp_ref_tot_1_1_1": "20000\n20100\n19850\n20050\n19980",
    },
    "S3 — Outlier detection demo (n=6)": {
        "gene_count": 1, "patient_count": 1, "num_ref_genes": 1,
        "ploidy": 2, "partition_vol": 0.85, "qc_min": 10000,
        "outlier_method": "Grubbs", "outlier_enabled": True,
        "description_tr": "1 hedef gen, n=6. Kontrol grubunda replikat 5 aykırı değer (kontaminasyon benzeri yüksek λ). Grubbs testiyle tespit edilir.",
        "description_en": "1 target gene, n=6. Replicate 5 in the control group is an outlier (contamination-like elevated λ). Detected by Grubbs' test.",
        "ctrl_tgt_pos_0": "1890\n1920\n1875\n1905\n5800\n1898",
        "ctrl_tgt_tot_0": "20000\n20100\n19850\n20050\n20000\n19980",
        "ctrl_ref_pos_0_0": "1895\n1915\n1870\n1900\n1888\n1892",
        "ctrl_ref_tot_0_0": "20000\n20100\n19850\n20050\n19990\n19980",
        "smp_tgt_pos_0_0": "3780\n3820\n3755\n3800\n3790",
        "smp_tgt_tot_0_0": "20000\n20100\n19850\n20050\n19980",
        "smp_ref_pos_0_0_0": "1895\n1915\n1870\n1900\n1892",
        "smp_ref_tot_0_0_0": "20000\n20100\n19850\n20050\n19980",
    },
    "S4 — Multi-group ANOVA (3 groups, n=5)": {
        "gene_count": 1, "patient_count": 3, "num_ref_genes": 1,
        "ploidy": 2, "partition_vol": 0.85, "qc_min": 10000,
        "outlier_method": "Grubbs", "outlier_enabled": True,
        "description_tr": "1 hedef gen, 3 hasta grubu, n=5. Grup 1: hafif kazanım, Grup 2: güçlü kazanım, Grup 3: değişim yok. Tek yönlü ANOVA + Tukey HSD.",
        "description_en": "1 target gene, 3 patient groups, n=5. Group 1: mild gain, Group 2: strong gain, Group 3: no change. One-way ANOVA + Tukey HSD.",
        "ctrl_tgt_pos_0": "1890\n1920\n1875\n1905\n1898",
        "ctrl_tgt_tot_0": "20000\n20100\n19850\n20050\n19980",
        "ctrl_ref_pos_0_0": "1895\n1915\n1870\n1900\n1892",
        "ctrl_ref_tot_0_0": "20000\n20100\n19850\n20050\n19980",
        # Group 1: mild gain (ratio ~1.3)
        "smp_tgt_pos_0_0": "2440\n2470\n2415\n2455\n2445",
        "smp_tgt_tot_0_0": "20000\n20100\n19850\n20050\n19980",
        "smp_ref_pos_0_0_0": "1895\n1915\n1870\n1900\n1892",
        "smp_ref_tot_0_0_0": "20000\n20100\n19850\n20050\n19980",
        # Group 2: strong gain (ratio ~2.5)
        "smp_tgt_pos_0_1": "4520\n4570\n4480\n4545\n4510",
        "smp_tgt_tot_0_1": "20000\n20100\n19850\n20050\n19980",
        "smp_ref_pos_0_1_0": "1895\n1915\n1870\n1900\n1892",
        "smp_ref_tot_0_1_0": "20000\n20100\n19850\n20050\n19980",
        # Group 3: no change
        "smp_tgt_pos_0_2": "1885\n1915\n1870\n1900\n1892",
        "smp_tgt_tot_0_2": "20000\n20100\n19850\n20050\n19980",
        "smp_ref_pos_0_2_0": "1895\n1915\n1870\n1900\n1892",
        "smp_ref_tot_0_2_0": "20000\n20100\n19850\n20050\n19980",
    },
    "S5 — QC & saturation demo (n=5)": {
        "gene_count": 1, "patient_count": 1, "num_ref_genes": 1,
        "ploidy": 2, "partition_vol": 0.85, "qc_min": 10000,
        "outlier_method": "Grubbs", "outlier_enabled": True,
        "description_tr": "1 hedef gen, n=5. Kontrol grubunda 1 replikat düşük partisyon (QC hatası), 1 replikat doygun (%100 pozitif). Kalite kontrol uyarılarını gösterir.",
        "description_en": "1 target gene, n=5. Control group has 1 low-partition replicate (QC failure) and 1 saturated replicate (100% positive). Demonstrates QC warnings.",
        "ctrl_tgt_pos_0": "1890\n1920\n500\n20000\n1898",
        "ctrl_tgt_tot_0": "20000\n20100\n3000\n20000\n19980",
        "ctrl_ref_pos_0_0": "1895\n1915\n505\n1901\n1892",
        "ctrl_ref_tot_0_0": "20000\n20100\n3000\n20050\n19980",
        "smp_tgt_pos_0_0": "3780\n3820\n3755\n3800\n3790",
        "smp_tgt_tot_0_0": "20000\n20100\n19850\n20050\n19980",
        "smp_ref_pos_0_0_0": "1895\n1915\n1870\n1900\n1892",
        "smp_ref_tot_0_0_0": "20000\n20100\n19850\n20050\n19980",
    },
}

# ═══════════════════════════════════════════════════════════════════════════════
# HEADER
# ═══════════════════════════════════════════════════════════════════════════════
_title_parts = _t['title'].split(' ', 1)
_header_emoji = _title_parts[0]
_header_title_text = _title_parts[1] if len(_title_parts) > 1 else _t['title']

st.markdown(
    f"""
    <div style="background:linear-gradient(90deg,#004d40,#00796b);
                color:white;padding:16px 18px;border-radius:8px;margin-bottom:8px;
                display:flex;align-items:center;gap:12px;">
        <span style="font-size:28px;line-height:1;flex-shrink:0;">{_header_emoji}</span>
        <div style="display:flex;flex-direction:column;justify-content:center;">
            <span style="font-size:20px;font-weight:800;line-height:1.3;">{_header_title_text}</span>
            <span style="font-size:11px;opacity:0.8;margin-top:3px;line-height:1.3;">{_t['subtitle']}</span>
        </div>
    </div>
    """,
    unsafe_allow_html=True
)

# ═══════════════════════════════════════════════════════════════════════════════
# SIDEBAR — User guide, links
# ═══════════════════════════════════════════════════════════════════════════════
# ═══════════════════════════════════════════════════════════════════════════════
# SIDEBAR — Example data loader
# ═══════════════════════════════════════════════════════════════════════════════
st.sidebar.markdown(
    f"<div style='font-size:12px;font-weight:600;color:#004d40;margin-bottom:2px;'>"
    f"{_t['sidebar_example_title']}</div>", unsafe_allow_html=True
)
selected_scenario = st.sidebar.selectbox(
    _t['sidebar_example_select'], options=["—"] + list(SCENARIOS.keys()),
    key="scenario_selector", label_visibility="collapsed"
)
if selected_scenario != "—":
    _sc = SCENARIOS[selected_scenario]
    st.sidebar.caption(_sc.get(f"description_{language_code}", _sc.get("description_en", "")))
    if st.sidebar.button(_t['sidebar_example_load_btn'], key="load_scenario_btn", use_container_width=True):
        for _key, _val in _sc.items():
            if _key.startswith("description"):
                continue
            st.session_state[_key] = _val
        st.sidebar.success(_t['sidebar_example_loaded'].format(s=selected_scenario))

# ═══════════════════════════════════════════════════════════════════════════════
# SIDEBAR — Instrument CSV import
# ═══════════════════════════════════════════════════════════════════════════════
with st.sidebar.expander(_t['csv_import_expander'], expanded=False):
    st.caption(_t['csv_import_description'])
    st.download_button(
        _t['download_example_csv'], data=EXAMPLE_CSV_MAIN_IMPORT.encode("utf-8"),
        file_name="example_import.csv", mime="text/csv", key="dl_example_main_import",
        use_container_width=True
    )
    csv_file = st.file_uploader(_t['csv_uploader'], type=["csv", "tsv", "txt"], key="csv_import_uploader")
    if csv_file is not None:
        _raw_bytes = csv_file.read()
        _raw_df, _detected_cols, _parse_err = parse_instrument_csv(_raw_bytes)
        if _parse_err:
            st.error(_t['csv_parse_error'].format(err=_parse_err))
        else:
            st.markdown(_t['csv_col_mapping_title'])
            _col_options = list(_raw_df.columns)
            _col_map = {}
            for _field, _label_key in [("sample", "csv_col_sample"), ("target", "csv_col_target"),
                                        ("positives", "csv_col_positives"), ("total", "csv_col_total")]:
                _detected = _detected_cols.get(_field)
                _status = _t['csv_col_auto_detected'] if _detected else _t['csv_col_manual_needed']
                _default_idx = _col_options.index(_detected) if _detected in _col_options else 0
                _col_map[_field] = st.selectbox(
                    f"{_t[_label_key]} — {_status}", options=_col_options,
                    index=_default_idx, key=f"csv_colmap_{_field}"
                )

            _std_df = build_standard_import_df(_raw_df, _col_map)
            if not _std_df.empty:
                with st.expander(_t['csv_preview_title'], expanded=False):
                    st.dataframe(_std_df, use_container_width=True)

                st.markdown(_t['csv_assay_assignment_title'])
                _unique_targets = sorted(_std_df["Target"].unique())
                _sel_target_assays = st.multiselect(
                    _t['csv_target_assays_label'], options=_unique_targets,
                    default=_unique_targets[:1] if _unique_targets else [],
                    key="csv_target_assays"
                )
                _remaining_targets = [t for t in _unique_targets if t not in _sel_target_assays]
                _sel_ref_assays = st.multiselect(
                    _t['csv_ref_assays_label'], options=_remaining_targets,
                    default=_remaining_targets[:1] if _remaining_targets else [],
                    key="csv_ref_assays"
                )

                st.markdown(_t['csv_group_assignment_title'])
                _unique_samples = sorted(_std_df["Sample"].unique())
                _ctrl_label_csv = st.text_input(
                    _t['csv_ctrl_label'],
                    value=_unique_samples[0] if _unique_samples else "",
                    key="csv_ctrl_label_input"
                )
                _n_pat_csv = st.number_input(
                    _t['csv_n_patient_groups'], min_value=1, max_value=10, value=1, step=1, key="csv_n_pat"
                )
                _patient_labels_csv = []
                for _pg in range(int(_n_pat_csv)):
                    _default_pat = _unique_samples[_pg + 1] if _pg + 1 < len(_unique_samples) else ""
                    _pat_lbl = st.text_input(
                        _t['csv_patient_label'].format(i=_pg + 1),
                        value=_default_pat, key=f"csv_pat_label_{_pg}"
                    )
                    _patient_labels_csv.append(_pat_lbl)

                if st.button(_t['csv_apply_btn'], key="csv_apply_btn", use_container_width=True):
                    if _sel_target_assays and _sel_ref_assays:
                        _n_filled = apply_csv_import_to_session(
                            _std_df, _sel_target_assays, _sel_ref_assays, _ctrl_label_csv, _patient_labels_csv
                        )
                        if _n_filled > 0:
                            st.success(_t['csv_apply_success'].format(n=_n_filled))
                        else:
                            st.warning(_t['csv_apply_warning'])
                    else:
                        st.warning(_t['csv_apply_warning'])

# ═══════════════════════════════════════════════════════════════════════════════
# SIDEBAR — Project save/load
# ═══════════════════════════════════════════════════════════════════════════════
with st.sidebar.expander(_t['project_expander'], expanded=False):
    st.caption(_t['project_description'])
    _project_data = export_project_state()
    st.download_button(
        _t['project_export_btn'],
        data=json.dumps(_project_data, ensure_ascii=False, indent=2).encode("utf-8"),
        file_name="absolutegene_project.json", mime="application/json",
        key="project_export_btn", use_container_width=True
    )
    project_file = st.file_uploader(_t['project_import_uploader'], type=["json"], key="project_import_uploader")
    if project_file is not None:
        try:
            _imported_project = json.loads(project_file.read().decode("utf-8"))
            _n_restored = import_project_state(_imported_project)
            st.success(_t['project_import_success'].format(n=_n_restored))
        except Exception as _e:
            st.error(_t['project_import_error'].format(err=str(_e)))

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

tab_data, tab_results, tab_batch, tab_vaf, tab_clinical, tab_crm, tab_report = st.tabs([
    f"📥 {_t['tab_data']}", f"📊 {_t['tab_results']}", f"🔬 {_t['tab_batch']}",
    f"🧬 {_t['tab_vaf']}", f"🩺 {_t['tab_clinical']}", f"🏭 {_t['tab_crm']}", f"📄 {_t['tab_report']}"
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
            num_ref_genes = st.number_input(_t['num_ref_genes'], min_value=1, max_value=10, step=1,
                                             key="num_ref_genes", help=_t['ref_gene_help'])
        sd_c4, sd_c5, sd_c6 = st.columns(3)
        with sd_c4:
            ploidy = st.number_input(_t['ploidy_label'], min_value=1, max_value=10,
                                      **({} if "ploidy" in st.session_state else {"value": 2}),
                                      step=1, key="ploidy", help=_t['ploidy_help'])
        with sd_c5:
            partition_vol_nl = st.number_input(
                _t['partition_vol_label'], min_value=0.01, max_value=100.0,
                **({} if "partition_vol" in st.session_state else {"value": 0.85}),
                step=0.01, format="%.2f", key="partition_vol", help=_t['partition_vol_help']
            )
        with sd_c6:
            qc_min_partitions = st.number_input(
                _t['qc_min_partitions'], min_value=100, max_value=100000,
                **({} if "qc_min" in st.session_state else {"value": 10000}),
                step=500, key="qc_min", help=_t['qc_min_partitions_help']
            )

    # ── Outlier Detection Settings ────────────────────────────────────────────
    with st.container(border=True):
        st.markdown(_t['outlier_section_title'])
        out_c1, out_c2 = st.columns([1, 2])
        with out_c1:
            outlier_enabled = st.checkbox(
                _t['outlier_enable'], key="outlier_enabled", help=_t['outlier_enable_help'],
                **({} if "outlier_enabled" in st.session_state else {"value": True})
            )
            outlier_method = st.radio(_t['outlier_method_label'], options=["Grubbs", "IQR"],
                                       key="outlier_method", help=_t['outlier_method_help'])
        with out_c2:
            if outlier_method == "Grubbs":
                grubbs_alpha = st.number_input(
                    _t['outlier_alpha_label'], min_value=0.01, max_value=0.10,
                    **({} if "grubbs_alpha" in st.session_state else {"value": 0.05}),
                    step=0.01, format="%.2f", key="grubbs_alpha"
                )
                iqr_multiplier = 1.5
            else:
                iqr_multiplier = st.number_input(
                    _t['outlier_iqr_label'], min_value=1.0, max_value=3.0,
                    **({} if "iqr_mult" in st.session_state else {"value": 1.5}),
                    step=0.25, format="%.2f", key="iqr_mult"
                )
                grubbs_alpha = 0.05

    st.divider()

    # ── Dilution / Dynamic Range Calculator ───────────────────────────────────
    with st.expander(_t['dilution_expander'], expanded=False):
        st.caption(_t['dilution_description'])
        dil_mode = st.radio(
            _t['dilution_mode_label'],
            options=[_t['dilution_mode_counts'], _t['dilution_mode_lambda']],
            key="dilution_mode", horizontal=True
        )
        dcol1, dcol2 = st.columns(2)
        if dil_mode == _t['dilution_mode_counts']:
            with dcol1:
                dil_pos = st.number_input(_t['dilution_positive_label'], min_value=0, value=1000, step=10, key="dil_pos")
            with dcol2:
                dil_tot = st.number_input(_t['dilution_total_label'], min_value=1, value=20000, step=100, key="dil_tot")
            _dil_lambda_input, _, _, _ = poisson_lambda(dil_pos, dil_tot)
        else:
            with dcol1:
                _dil_lambda_input = st.number_input(_t['dilution_lambda_label'], min_value=0.001, value=1.6,
                                                     step=0.01, format="%.3f", key="dil_lambda_direct")

        dil_target = st.number_input(_t['dilution_target_label'], min_value=0.1, max_value=3.0,
                                      value=1.6, step=0.1, format="%.2f", key="dil_target_lambda")

        if st.button(_t['dilution_calc_btn'], key="dil_calc_btn"):
            _dil_result = recommend_dilution(_dil_lambda_input, dil_target)
            if _dil_result is None:
                st.error("Invalid input.")
            else:
                _lam = _dil_result["current_lambda"]
                _factor = _dil_result["dilution_factor"]
                if _dil_result["status"] == "optimal":
                    st.success(_t['dilution_result_optimal'].format(lam=_lam))
                elif _dil_result["status"] == "too_low":
                    st.warning(_t['dilution_result_too_low'].format(lam=_lam, factor=_factor))
                elif _dil_result["status"] == "too_high":
                    st.warning(_t['dilution_result_too_high'].format(lam=_lam, factor=_factor))
                else:
                    st.error(_t['dilution_result_saturated'].format(lam=_lam, factor=_factor))
                st.caption(_t['dilution_factor_note'])

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

    def _ta(label, key):
        """text_area wrapper that avoids the Streamlit session-state/value conflict warning."""
        kwargs = {} if key in st.session_state else {"value": ""}
        return st.text_area(label, key=key, **kwargs)

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
            ctrl_target_pos_txt = _ta(
                f"Control {i+1} — {_t['positive_partitions']} ({_t['target_gene']})", f"ctrl_tgt_pos_{i}"
            )
        with cc2:
            ctrl_target_tot_txt = _ta(
                f"Control {i+1} — {_t['total_partitions']} ({_t['target_gene']})", f"ctrl_tgt_tot_{i}"
            )

        ctrl_ref_pos_txts, ctrl_ref_tot_txts = [], []
        for r in range(num_ref_genes):
            rc1, rc2 = st.columns(2)
            ref_lbl = f"{_t['reference_gene']} {r+1}" if num_ref_genes > 1 else _t['reference_gene']
            with rc1:
                rp = _ta(f"Control {i+1} — {_t['positive_partitions']} ({ref_lbl})", f"ctrl_ref_pos_{i}_{r}")
            with rc2:
                rt = _ta(f"Control {i+1} — {_t['total_partitions']} ({ref_lbl})", f"ctrl_ref_tot_{i}_{r}")
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
                "__conc__": round(row["target_lambda"] / partition_vol_nl * 1000.0, 2) if not np.isnan(row["target_lambda"]) else None,
                "Outlier Excluded": status_label,
                "__used__": row["used"],
            })

        avg_ctrl_ratio = float(np.mean(ctrl_result["ratio_kept"])) if len(ctrl_result["ratio_kept"]) > 0 else None

        # ── NTC / LOD-LOQ (optional) ──────────────────────────────────────────
        lod_loq_result = None
        with st.expander(_t['ntc_expander'], expanded=False):
            st.caption(_t['ntc_description'])
            ntc_c1, ntc_c2 = st.columns(2)
            with ntc_c1:
                ntc_pos_txt = _ta(f"{_t['ntc_positive_label']} — {_t['target_gene']} {i+1}", f"ntc_pos_{i}")
            with ntc_c2:
                ntc_tot_txt = _ta(f"{_t['ntc_total_label']} — {_t['target_gene']} {i+1}", f"ntc_tot_{i}")

            ntc_pos_arr = parse_input_data(ntc_pos_txt)
            ntc_tot_arr = parse_input_data(ntc_tot_txt)
            if len(ntc_pos_arr) > 0 and len(ntc_tot_arr) > 0:
                lod_loq_result = compute_lod_loq(ntc_pos_arr, ntc_tot_arr, partition_vol_nl)
                if lod_loq_result is not None:
                    _gene_label_lod = f"{_t['target_gene']} {i+1}"
                    st.markdown(f"**{_t['lod_result_title'].format(gene=_gene_label_lod)}**")
                    lod_c1, lod_c2 = st.columns(2)
                    lod_c1.metric(_t['lod_label'], f"{lod_loq_result['lod_conc']:.2f}")
                    lod_c2.metric(_t['loq_label'], f"{lod_loq_result['loq_conc']:.2f}")
                    if lod_loq_result["contamination"]:
                        st.warning(_t['ntc_contamination_warning'].format(
                            lam=lod_loq_result["ntc_lambda"],
                            pos=int(lod_loq_result["ntc_pos_total"]), tot=int(lod_loq_result["ntc_tot_total"])
                        ))
                    else:
                        st.info(_t['ntc_zero_note'].format(n=int(lod_loq_result["ntc_tot_total"])))
                    st.caption(_t['loq_heuristic_note'])

        # ── Patient groups ────────────────────────────────────────────────────
        for j in range(num_patient_groups):
            st.markdown(f"**{_t['patient_group']} {j+1} — {_t['target_gene']} {i+1}**")
            pc1, pc2 = st.columns(2)
            with pc1:
                smp_target_pos_txt = _ta(
                    f"Group {j+1} — {_t['positive_partitions']} ({_t['target_gene']} {i+1})", f"smp_tgt_pos_{i}_{j}"
                )
            with pc2:
                smp_target_tot_txt = _ta(
                    f"Group {j+1} — {_t['total_partitions']} ({_t['target_gene']} {i+1})", f"smp_tgt_tot_{i}_{j}"
                )

            smp_ref_pos_txts, smp_ref_tot_txts = [], []
            for r in range(num_ref_genes):
                rc1, rc2 = st.columns(2)
                ref_lbl = f"{_t['reference_gene']} {r+1}" if num_ref_genes > 1 else _t['reference_gene']
                with rc1:
                    rp = _ta(f"Group {j+1} — {_t['positive_partitions']} ({ref_lbl})", f"smp_ref_pos_{i}_{j}_{r}")
                with rc2:
                    rt = _ta(f"Group {j+1} — {_t['total_partitions']} ({ref_lbl})", f"smp_ref_tot_{i}_{j}_{r}")
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
                    "__conc__": round(row["target_lambda"] / partition_vol_nl * 1000.0, 2) if not np.isnan(row["target_lambda"]) else None,
                    "Outlier Excluded": status_label,
                    "__used__": row["used"],
                })

            avg_smp_ratio = float(np.mean(smp_result["ratio_kept"]))
            fold_change = avg_smp_ratio / avg_ctrl_ratio if avg_ctrl_ratio else float("nan")
            cn_ctrl = ploidy * avg_ctrl_ratio if avg_ctrl_ratio else float("nan")
            cn_smp = ploidy * avg_smp_ratio

            # ── Concentration (copies/µL) with 95% CI (replicate-based) ─────────
            conc_ctrl_arr = lambda_to_conc(ctrl_result["lam_target_kept"], partition_vol_nl)
            conc_smp_arr = lambda_to_conc(smp_result["lam_target_kept"], partition_vol_nl)
            conc_ctrl_mean, conc_ctrl_lo, conc_ctrl_hi, conc_ctrl_cv, _ = mean_ci(conc_ctrl_arr)
            conc_smp_mean, conc_smp_lo, conc_smp_hi, conc_smp_cv, _ = mean_ci(conc_smp_arr)

            # ── LOD/LOQ QC flag (if NTC data provided for this gene) ─────────────
            lod_loq_flag = "—"
            if lod_loq_result is not None:
                if conc_smp_mean < lod_loq_result["lod_conc"]:
                    lod_loq_flag = _t['below_lod_flag']
                elif conc_smp_mean < lod_loq_result["loq_conc"]:
                    lod_loq_flag = _t['between_lod_loq_flag']
                else:
                    lod_loq_flag = _t['above_loq_flag']

            if fold_change >= 1.5:
                regulation = _t['upregulated']
            elif fold_change <= 0.67:
                regulation = _t['downregulated']
            else:
                regulation = _t['no_change']

            st.markdown(f"#### {_t['method_comparison']} — {_t['target_gene']} {i+1} / {_t['patient_group']} {j+1}")
            rcol1, rcol2, rcol3, rcol4 = st.columns(4)
            rcol1.metric(_t['ratio_col'], f"{avg_smp_ratio:.4f}")
            rcol2.metric(_t['cn_col'], f"{cn_smp:.3f}")
            rcol3.metric(_t['fc_col'], f"{fold_change:.4f}", delta=regulation)
            _conc_ci_txt = (f"95% CI: {conc_smp_lo:.1f}–{conc_smp_hi:.1f}"
                             if not np.isnan(conc_smp_lo) else "n<2, no CI")
            rcol4.metric(_t['conc_col'], f"{conc_smp_mean:.1f}", delta=_conc_ci_txt, delta_color="off")
            if lod_loq_result is not None:
                st.caption(f"{_t['lod_qc_col']}: {lod_loq_flag}")

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
                "__conc_ctrl__": conc_ctrl_mean, "__conc_ctrl_lo__": conc_ctrl_lo, "__conc_ctrl_hi__": conc_ctrl_hi,
                "__conc_smp__": conc_smp_mean, "__conc_smp_lo__": conc_smp_lo, "__conc_smp_hi__": conc_smp_hi,
                "__conc_smp_cv__": conc_smp_cv, "__lod_loq_flag__": lod_loq_flag,
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

    # ── QC Summary Panel ──────────────────────────────────────────────────────
    if input_values_table:
        _n_total_rep = len(input_values_table)
        _n_qc_fail = sum(1 for r in input_values_table if r["Outlier Excluded"] == _t['qc_fail'])
        _n_saturated = sum(1 for r in input_values_table if r["Outlier Excluded"] == _t['qc_saturated'])
        _n_outlier = sum(1 for r in input_values_table if str(r["Outlier Excluded"]).startswith(_t['outlier_excluded_yes']))
        _n_below_lod = sum(1 for r in data if r.get("__lod_loq_flag__") in (_t['below_lod_flag'], _t['between_lod_loq_flag']))
        _n_high_cv = sum(1 for r in data if r.get("__conc_smp_cv__") is not None
                          and not np.isnan(r["__conc_smp_cv__"]) and r["__conc_smp_cv__"] > 25)

        with st.container(border=True):
            st.markdown(f"#### {_t['qc_panel_title']}")
            qc_c1, qc_c2, qc_c3, qc_c4, qc_c5, qc_c6 = st.columns(6)
            qc_c1.metric(_t['qc_panel_total'], _n_total_rep)
            qc_c2.metric(_t['qc_panel_qc_fail'], _n_qc_fail)
            qc_c3.metric(_t['qc_panel_saturated'], _n_saturated)
            qc_c4.metric(_t['qc_panel_outlier'], _n_outlier)
            qc_c5.metric(_t['qc_panel_below_lod'], _n_below_lod)
            qc_c6.metric(_t['qc_panel_high_cv'], _n_high_cv)

            _n_flags_total = _n_qc_fail + _n_saturated + _n_outlier + _n_below_lod + _n_high_cv
            _flag_ratio = _n_flags_total / max(_n_total_rep, 1)
            if _flag_ratio == 0:
                st.success(_t['qc_panel_verdict_good'])
            elif _flag_ratio < 0.15:
                st.warning(_t['qc_panel_verdict_caution'])
            else:
                st.error(_t['qc_panel_verdict_poor'])
        st.markdown("---")
    else:
        st.info(_t['qc_panel_no_data'])

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
            "__lambda__": _t['lambda_col'], "__conc__": _t['conc_col'], "Outlier Excluded": "Status",
        }
        display_df = pd.DataFrame(input_values_table).drop(columns=["__used__"]).rename(columns=_rename)
        st.dataframe(display_df, use_container_width=True)
        csv = display_df.to_csv(index=False).encode("utf-8")
        st.download_button(_t['download_csv'], data=csv, file_name="dpcr_input_data.csv", mime="text/csv", key="dl_input_csv")

    # ── Results summary table ─────────────────────────────────────────────────
    if data:
        st.subheader(_t['nil_mine'])
        _res_display_rows = []
        for r in data:
            _res_display_rows.append({
                _t['target_gene']: r["__gene__"], "Group": r["__group__"],
                f"{_t['ratio_col']} (Control)": r["__ratio_ctrl__"], f"{_t['ratio_col']} (Sample)": r["__ratio_smp__"],
                f"{_t['cn_col']} (Control)": r["__cn_ctrl__"], f"{_t['cn_col']} (Sample)": r["__cn_smp__"],
                f"{_t['conc_col']} (Control)": (f"{r['__conc_ctrl__']:.1f}" if not np.isnan(r['__conc_ctrl__']) else "—"),
                f"{_t['conc_col']} (Sample)": (f"{r['__conc_smp__']:.1f}" if not np.isnan(r['__conc_smp__']) else "—"),
                "95% CI (Sample, copies/µL)": (
                    f"{r['__conc_smp_lo__']:.1f}–{r['__conc_smp_hi__']:.1f}"
                    if not np.isnan(r['__conc_smp_lo__']) else "n<2"
                ),
                "CV% (Sample)": (f"{r['__conc_smp_cv__']:.1f}" if not np.isnan(r['__conc_smp_cv__']) else "—"),
                _t['lod_qc_col']: r.get("__lod_loq_flag__", "—"),
                _t['fc_col']: r["__fc__"], _t['regulation_status']: r["__regulation__"],
                "n (Control)": r["__n_ctrl__"], "n (Sample)": r["__n_smp__"],
            })
        res_df = pd.DataFrame(_res_display_rows)
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
def create_simple_pdf(report_title, subtitle, description, summary_rows, table_header,
                       table_rows, chart_png_bytes, chart_caption, footer_note, references=None):
    """
    Generic lightweight PDF builder shared by the Batch Screening and VAF
    Calculator tabs. Produces: title/subtitle, description paragraph,
    a summary key-value table, a results table, an optional chart image,
    a footer note, and an optional references list.
    """
    buffer = BytesIO()
    doc = SimpleDocTemplate(buffer, pagesize=letter, leftMargin=50, rightMargin=50, topMargin=60, bottomMargin=50)
    elements = []
    styles = getSampleStyleSheet()

    title_style = ParagraphStyle('RT', parent=styles['Title'], fontName=PDF_FONT_BOLD, fontSize=18,
                                  textColor=colors.HexColor('#004d40'), spaceAfter=6, alignment=1)
    sub_style = ParagraphStyle('RS', parent=styles['Normal'], fontName=PDF_FONT, fontSize=10,
                                textColor=colors.HexColor('#555555'), spaceAfter=4, alignment=1)
    body_style = ParagraphStyle('BD', parent=styles['Normal'], fontName=PDF_FONT, fontSize=9, leading=13, spaceAfter=8)
    small_style = ParagraphStyle('SM', parent=styles['Normal'], fontName=PDF_FONT, fontSize=8, leading=11,
                                  textColor=colors.HexColor('#444444'))
    caption_style = ParagraphStyle('CA', parent=styles['Normal'], fontName=PDF_FONT, fontSize=8,
                                    textColor=colors.HexColor('#666666'), alignment=1, spaceAfter=6)

    def hr():
        return HRFlowable(width="100%", thickness=0.5, color=colors.HexColor('#cccccc'), spaceAfter=8, spaceBefore=4)

    def make_table(rows, col_widths=None, header=True):
        if not rows:
            return Spacer(1, 1)
        styled_rows = []
        for ri, row in enumerate(rows):
            styled_row = []
            for cell in row:
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

    elements.append(Spacer(1, 30))
    elements.append(Paragraph(safe_str(report_title), title_style))
    elements.append(Paragraph(safe_str(subtitle), sub_style))
    elements.append(Paragraph(safe_str(f"Generated: {now}"), sub_style))
    elements.append(Spacer(1, 14))
    elements.append(hr())
    elements.append(Paragraph(safe_str(description), body_style))

    if summary_rows:
        elements.append(make_table(summary_rows, col_widths=[260, 200]))
        elements.append(Spacer(1, 14))

    if table_header and table_rows:
        cw = (letter[0] - 100) / len(table_header)
        elements.append(make_table([table_header] + table_rows, col_widths=[cw] * len(table_header)))
        elements.append(Spacer(1, 10))

    if chart_png_bytes:
        elements.append(RLImage(BytesIO(chart_png_bytes), width=460, height=250))
        if chart_caption:
            elements.append(Paragraph(safe_str(chart_caption), caption_style))

    if references:
        elements.append(Spacer(1, 12))
        elements.append(hr())
        for ref in references:
            elements.append(Paragraph(safe_str(f"• {ref}"), small_style))
            elements.append(Spacer(1, 3))

    elements.append(Spacer(1, 14))
    elements.append(hr())
    elements.append(Paragraph(safe_str(footer_note), small_style))

    doc.build(elements)
    buffer.seek(0)
    return buffer


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
    res_cols = ["Gene", "Group", "Ratio (Ctrl)", "Ratio (Sample)", "CN (Ctrl)", "CN (Sample)",
                "Conc. Sample (copies/µL)", "95% CI (copies/µL)", "Fold Change", "Regulation", "LOD/LOQ Status"]
    res_rows = [res_cols]
    for r in results:
        conc_ci = (f"{r['__conc_smp_lo__']:.1f}\u2013{r['__conc_smp_hi__']:.1f}"
                   if r.get('__conc_smp_lo__') is not None and not np.isnan(r['__conc_smp_lo__']) else "n<2")
        res_rows.append([
            r["__gene__"], r["__group__"],
            f"{r['__ratio_ctrl__']:.4f}" if r['__ratio_ctrl__'] is not None else "—",
            f"{r['__ratio_smp__']:.4f}",
            f"{r['__cn_ctrl__']:.3f}" if r['__cn_ctrl__'] is not None else "—",
            f"{r['__cn_smp__']:.3f}",
            f"{r.get('__conc_smp__', float('nan')):.1f}" if r.get('__conc_smp__') is not None and not np.isnan(r.get('__conc_smp__', np.nan)) else "—",
            conc_ci,
            f"{r['__fc__']:.4f}", r["__regulation__"],
            r.get("__lod_loq_flag__", "—"),
        ])
    cw11r = (letter[0] - 100) / 11
    elements.append(make_table(res_rows, col_widths=[cw11r] * 11))
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


# ═══════════════════════════════════════════════════════════════════════════════
# BATCH SCREENING TAB
# ═══════════════════════════════════════════════════════════════════════════════
with tab_batch:
    st.markdown(f"### {_t['batch_title']}")
    st.caption(_t['batch_description'])
    st.download_button(
        _t['download_example_csv'], data=EXAMPLE_CSV_BATCH_SCREENING.encode("utf-8"),
        file_name="example_batch_screening.csv", mime="text/csv", key="dl_example_batch"
    )
    st.markdown("---")

    batch_file = st.file_uploader(_t['batch_uploader'], type=["csv", "tsv", "txt"], key="batch_uploader")

    if batch_file is not None:
        _batch_raw_bytes = batch_file.read()
        _batch_raw_df, _batch_detected, _batch_err = parse_instrument_csv(_batch_raw_bytes)
        if _batch_err:
            st.error(_t['csv_parse_error'].format(err=_batch_err))
        else:
            _batch_col_options = list(_batch_raw_df.columns)
            _batch_col_map = {}
            bcol1, bcol2, bcol3, bcol4 = st.columns(4)
            for _field, _label_key, _col_widget in [
                ("sample", "csv_col_sample", bcol1), ("target", "csv_col_target", bcol2),
                ("positives", "csv_col_positives", bcol3), ("total", "csv_col_total", bcol4)
            ]:
                _detected = _batch_detected.get(_field)
                _default_idx = _batch_col_options.index(_detected) if _detected in _batch_col_options else 0
                with _col_widget:
                    _batch_col_map[_field] = st.selectbox(
                        _t[_label_key], options=_batch_col_options,
                        index=_default_idx, key=f"batch_colmap_{_field}"
                    )

            _batch_std_df = build_standard_import_df(_batch_raw_df, _batch_col_map)

            if not _batch_std_df.empty:
                st.markdown(_t['batch_assay_title'])
                _batch_unique_targets = sorted(_batch_std_df["Target"].unique())
                bacol1, bacol2 = st.columns(2)
                with bacol1:
                    _batch_target_assay = st.selectbox(
                        _t['batch_target_label'], options=_batch_unique_targets, key="batch_target_assay"
                    )
                with bacol2:
                    _batch_ref_options = [t for t in _batch_unique_targets if t != _batch_target_assay]
                    _batch_ref_assays = st.multiselect(
                        _t['batch_ref_label'], options=_batch_ref_options,
                        default=_batch_ref_options[:1] if _batch_ref_options else [],
                        key="batch_ref_assays"
                    )

                st.markdown(_t['batch_settings_title'])
                bscol1, bscol2 = st.columns(2)
                with bscol1:
                    _batch_expected_ratio = st.number_input(
                        _t['batch_expected_ratio'], min_value=0.01, max_value=100.0,
                        value=1.0, step=0.1, key="batch_expected_ratio"
                    )
                with bscol2:
                    _batch_ploidy = st.number_input(
                        _t['batch_ploidy'], min_value=1, max_value=10, value=2, step=1, key="batch_ploidy"
                    )

                if st.button(_t['batch_run_btn'], key="batch_run_btn", use_container_width=True) and _batch_ref_assays:
                    _batch_results = pool_and_compute_batch(
                        _batch_std_df, _batch_target_assay, _batch_ref_assays, partition_vol_nl
                    )
                    st.session_state["_batch_results_cache"] = _batch_results
                    st.session_state["_batch_expected_ratio_cache"] = _batch_expected_ratio
                    st.session_state["_batch_ploidy_cache"] = _batch_ploidy
                    st.session_state["_batch_target_assay_cache"] = _batch_target_assay

    # ── Display cached results (persists across reruns/tab switches) ───────────
    _batch_results = st.session_state.get("_batch_results_cache")
    if _batch_results:
        _exp_ratio = st.session_state.get("_batch_expected_ratio_cache", 1.0)
        _b_ploidy = st.session_state.get("_batch_ploidy_cache", 2)
        _b_target_assay = st.session_state.get("_batch_target_assay_cache", "—")

        st.markdown(f"#### {_t['batch_results_title']}")
        _batch_rows = []
        _n_flagged = 0
        for r in _batch_results:
            if np.isnan(r["ci_low"]) or np.isnan(r["ci_high"]):
                classification = "—"
            elif _exp_ratio < r["ci_low"]:
                classification = _t['batch_class_gain']
                _n_flagged += 1
            elif _exp_ratio > r["ci_high"]:
                classification = _t['batch_class_loss']
                _n_flagged += 1
            else:
                classification = _t['batch_class_normal']

            _batch_rows.append({
                _t['batch_col_sample']: r["Sample"],
                _t['batch_col_lambda_t']: round(r["lambda_t"], 5),
                _t['batch_col_lambda_r']: round(r["lambda_r"], 5),
                _t['batch_col_ratio']: round(r["ratio"], 4),
                _t['batch_col_ci']: (f"{r['ci_low']:.4f}\u2013{r['ci_high']:.4f}"
                                       if not np.isnan(r["ci_low"]) else "—"),
                _t['batch_col_cn']: round(_b_ploidy * r["ratio"], 3),
                _t['batch_col_conc']: round(r["conc"], 1),
                _t['batch_col_class']: classification,
            })

        st.info(_t['batch_n_flagged'].format(n=_n_flagged, total=len(_batch_results)))
        batch_df = pd.DataFrame(_batch_rows)
        st.dataframe(batch_df, use_container_width=True)

        batch_csv = batch_df.to_csv(index=False).encode("utf-8")
        st.download_button(_t['batch_download_btn'], data=batch_csv,
                            file_name="batch_screening_results.csv", mime="text/csv", key="batch_dl_csv")

        # ── Screening plot ────────────────────────────────────────────────────
        fig_batch = go.Figure()
        _colors_map = {
            _t['batch_class_gain']: "#e53935", _t['batch_class_loss']: "#1e88e5",
            _t['batch_class_normal']: "#43a047", "—": "#9e9e9e"
        }
        _sample_names = [r["Sample"] for r in _batch_results]
        _ratios = [r["ratio"] for r in _batch_results]
        _err_low = [max(r["ratio"] - r["ci_low"], 0) if not np.isnan(r["ci_low"]) else 0 for r in _batch_results]
        _err_high = [r["ci_high"] - r["ratio"] if not np.isnan(r["ci_high"]) else 0 for r in _batch_results]
        _point_colors = [_colors_map.get(row[_t['batch_col_class']], "#9e9e9e") for row in _batch_rows]

        fig_batch.add_trace(go.Scatter(
            x=_sample_names, y=_ratios, mode="markers",
            marker=dict(color=_point_colors, size=9),
            error_y=dict(type="data", symmetric=False, array=_err_high, arrayminus=_err_low),
            name="Ratio (95% CI)"
        ))
        fig_batch.add_hline(y=_exp_ratio, line_dash="dash", line_color="black",
                            annotation_text=f"Expected = {_exp_ratio}")
        fig_batch.update_layout(
            title=_t['batch_chart_title'], xaxis_title=_t['batch_col_sample'],
            yaxis_title=_t['batch_col_ratio'], height=420, showlegend=False
        )
        st.plotly_chart(fig_batch, use_container_width=True, key="batch_screening_chart")

        # ── PDF report ────────────────────────────────────────────────────────
        if st.button(_t['batch_pdf_btn'], key="batch_pdf_btn"):
            try:
                fig_mpl, ax_mpl = plt.subplots(figsize=(7, 3.5))
                _mpl_colors = {
                    _t['batch_class_gain']: "#e53935", _t['batch_class_loss']: "#1e88e5",
                    _t['batch_class_normal']: "#43a047", "—": "#9e9e9e"
                }
                _mpl_point_colors = [_mpl_colors.get(row[_t['batch_col_class']], "#9e9e9e") for row in _batch_rows]
                ax_mpl.errorbar(range(len(_sample_names)), _ratios, yerr=[_err_low, _err_high],
                                 fmt='none', ecolor='gray', capsize=3, zorder=2)
                ax_mpl.scatter(range(len(_sample_names)), _ratios, c=_mpl_point_colors, s=50, zorder=3)
                ax_mpl.axhline(y=_exp_ratio, color='black', linestyle='--', linewidth=1)
                ax_mpl.set_xticks(range(len(_sample_names)))
                ax_mpl.set_xticklabels(_sample_names, rotation=45, ha='right', fontsize=7)
                ax_mpl.set_ylabel(_t['batch_col_ratio'], fontsize=9)
                ax_mpl.set_title(_t['batch_chart_title'], fontsize=10, fontweight='bold')
                ax_mpl.spines['top'].set_visible(False)
                ax_mpl.spines['right'].set_visible(False)
                plt.tight_layout()
                _img_buf = BytesIO()
                plt.savefig(_img_buf, format='png', dpi=150, bbox_inches='tight')
                plt.close()
                _chart_bytes = _img_buf.getvalue()
            except Exception:
                _chart_bytes = None

            _batch_summary_rows = [
                ["Parameter", "Value"],
                ["Target assay", _b_target_assay],
                ["Samples screened", str(len(_batch_results))],
                ["Expected ratio", f"{_exp_ratio}"],
                ["Reference ploidy", str(_b_ploidy)],
                ["Samples flagged (significant deviation)", str(_n_flagged)],
            ]
            _batch_pdf_table_header = [_t['batch_col_sample'], _t['batch_col_lambda_t'], _t['batch_col_lambda_r'],
                                        _t['batch_col_ratio'], _t['batch_col_ci'], _t['batch_col_cn'],
                                        _t['batch_col_conc'], _t['batch_col_class']]
            _batch_pdf_table_rows = [[str(v) for v in row.values()] for row in _batch_rows]

            batch_pdf_buffer = create_simple_pdf(
                report_title="AbsoluteGene",
                subtitle=_t['batch_pdf_report'],
                description=_t['batch_pdf_description'],
                summary_rows=_batch_summary_rows,
                table_header=_batch_pdf_table_header,
                table_rows=_batch_pdf_table_rows,
                chart_png_bytes=_chart_bytes,
                chart_caption="Figure. Ratio (95% CI) per sample. Dashed line = expected (no-change) ratio.",
                footer_note="AbsoluteGene — For research and educational use only. Not validated for clinical diagnostic purposes.",
                references=[
                    "Poisson-based dPCR quantification: Dube S, Qin J, Ramakrishnan R (2008). PLoS ONE, 3(8), e2876.",
                    "digital MIQE guidelines: Huggett JF et al. (2013). Clin Chem, 59(6), 892-902.",
                ]
            )
            st.download_button(
                "⬇️ " + _t['batch_pdf_report'], data=batch_pdf_buffer,
                file_name="batch_screening_report.pdf", mime="application/pdf", key="batch_pdf_dl"
            )
    else:
        st.info(_t['batch_no_data'])


# ═══════════════════════════════════════════════════════════════════════════════
# VAF / MUTATION FRACTION CALCULATOR TAB
# ═══════════════════════════════════════════════════════════════════════════════
with tab_vaf:
    st.markdown(f"### {_t['vaf_title']}")
    st.caption(_t['vaf_description'])
    st.info(_t['vaf_method_note'])
    st.download_button(
        _t['download_example_csv'], data=EXAMPLE_CSV_VAF.encode("utf-8"),
        file_name="example_vaf.csv", mime="text/csv", key="dl_example_vaf"
    )
    st.markdown("---")

    vaf_file = st.file_uploader(_t['vaf_uploader'], type=["csv", "tsv", "txt"], key="vaf_uploader")

    if vaf_file is not None:
        _vaf_raw_bytes = vaf_file.read()
        _vaf_raw_df, _vaf_detected, _vaf_err = parse_instrument_csv(_vaf_raw_bytes)
        if _vaf_err:
            st.error(_t['csv_parse_error'].format(err=_vaf_err))
        else:
            _vaf_col_options = list(_vaf_raw_df.columns)
            _vaf_col_map = {}
            vcol1, vcol2, vcol3, vcol4 = st.columns(4)
            for _field, _label_key, _col_widget in [
                ("sample", "csv_col_sample", vcol1), ("target", "csv_col_target", vcol2),
                ("positives", "csv_col_positives", vcol3), ("total", "csv_col_total", vcol4)
            ]:
                _detected = _vaf_detected.get(_field)
                _default_idx = _vaf_col_options.index(_detected) if _detected in _vaf_col_options else 0
                with _col_widget:
                    _vaf_col_map[_field] = st.selectbox(
                        _t[_label_key], options=_vaf_col_options,
                        index=_default_idx, key=f"vaf_colmap_{_field}"
                    )

            _vaf_std_df = build_standard_import_df(_vaf_raw_df, _vaf_col_map)

            if not _vaf_std_df.empty:
                st.markdown(_t['vaf_assay_title'])
                _vaf_unique_targets = sorted(_vaf_std_df["Target"].unique())
                vacol1, vacol2 = st.columns(2)
                with vacol1:
                    _vaf_mutant_assay = st.selectbox(
                        _t['vaf_mutant_label'], options=_vaf_unique_targets, key="vaf_mutant_assay"
                    )
                with vacol2:
                    _vaf_wt_options = [t for t in _vaf_unique_targets if t != _vaf_mutant_assay]
                    _vaf_wt_assay = st.selectbox(
                        _t['vaf_wt_label'], options=_vaf_wt_options, key="vaf_wt_assay"
                    ) if _vaf_wt_options else None

                if st.button(_t['vaf_run_btn'], key="vaf_run_btn", use_container_width=True) and _vaf_wt_assay:
                    _vaf_results = pool_and_compute_vaf(
                        _vaf_std_df, _vaf_mutant_assay, _vaf_wt_assay, partition_vol_nl
                    )
                    st.session_state["_vaf_results_cache"] = _vaf_results
                    st.session_state["_vaf_mutant_assay_cache"] = _vaf_mutant_assay
                    st.session_state["_vaf_wt_assay_cache"] = _vaf_wt_assay

    # ── Display cached results ──────────────────────────────────────────────────
    _vaf_results = st.session_state.get("_vaf_results_cache")
    if _vaf_results:
        st.markdown(f"#### {_t['vaf_results_title']}")
        _vaf_rows = []
        _n_detected = 0
        for r in _vaf_results:
            if r["detected"]:
                _n_detected += 1
            _vaf_rows.append({
                _t['vaf_col_sample']: r["Sample"],
                _t['vaf_col_lambda_mut']: round(r["lambda_mut"], 5),
                _t['vaf_col_lambda_wt']: round(r["lambda_wt"], 5),
                _t['vaf_col_fa']: round(r["fa"] * 100, 3),
                _t['vaf_col_ci']: (f"{r['ci_low']*100:.3f}\u2013{r['ci_high']*100:.3f}"
                                     if not np.isnan(r["ci_low"]) else "—"),
                _t['vaf_col_conc_mut']: round(r["conc_mut"], 2),
                _t['vaf_col_detected']: _t['vaf_detected_yes'] if r["detected"] else _t['vaf_detected_no'],
            })

        st.info(_t['vaf_n_detected'].format(n=_n_detected, total=len(_vaf_results)))
        vaf_df = pd.DataFrame(_vaf_rows)
        st.dataframe(vaf_df, use_container_width=True)

        vaf_csv = vaf_df.to_csv(index=False).encode("utf-8")
        st.download_button(_t['vaf_download_btn'], data=vaf_csv,
                            file_name="vaf_results.csv", mime="text/csv", key="vaf_dl_csv")

        # ── VAF plot ───────────────────────────────────────────────────────────
        fig_vaf = go.Figure()
        _vaf_samples = [r["Sample"] for r in _vaf_results]
        _vaf_pct = [r["fa"] * 100 for r in _vaf_results]
        _vaf_err_low = [max((r["fa"] - r["ci_low"]) * 100, 0) if not np.isnan(r["ci_low"]) else 0 for r in _vaf_results]
        _vaf_err_high = [(r["ci_high"] - r["fa"]) * 100 if not np.isnan(r["ci_high"]) else 0 for r in _vaf_results]
        _vaf_colors = ["#e53935" if r["detected"] else "#9e9e9e" for r in _vaf_results]

        fig_vaf.add_trace(go.Bar(
            x=_vaf_samples, y=_vaf_pct, marker_color=_vaf_colors,
            error_y=dict(type="data", symmetric=False, array=_vaf_err_high, arrayminus=_vaf_err_low),
            name="FA%"
        ))
        fig_vaf.update_layout(
            title=_t['vaf_chart_title'], xaxis_title=_t['vaf_col_sample'],
            yaxis_title=_t['vaf_col_fa'], height=420, showlegend=False
        )
        st.plotly_chart(fig_vaf, use_container_width=True, key="vaf_chart")

        # ── PDF report ────────────────────────────────────────────────────────
        if st.button(_t['vaf_pdf_btn'], key="vaf_pdf_btn"):
            try:
                fig_mpl, ax_mpl = plt.subplots(figsize=(7, 3.5))
                _mpl_vaf_colors = ["#e53935" if r["detected"] else "#9e9e9e" for r in _vaf_results]
                ax_mpl.bar(range(len(_vaf_samples)), _vaf_pct, color=_mpl_vaf_colors,
                           yerr=[_vaf_err_low, _vaf_err_high], capsize=3)
                ax_mpl.set_xticks(range(len(_vaf_samples)))
                ax_mpl.set_xticklabels(_vaf_samples, rotation=45, ha='right', fontsize=7)
                ax_mpl.set_ylabel(_t['vaf_col_fa'], fontsize=9)
                ax_mpl.set_title(_t['vaf_chart_title'], fontsize=10, fontweight='bold')
                ax_mpl.spines['top'].set_visible(False)
                ax_mpl.spines['right'].set_visible(False)
                plt.tight_layout()
                _img_buf = BytesIO()
                plt.savefig(_img_buf, format='png', dpi=150, bbox_inches='tight')
                plt.close()
                _vaf_chart_bytes = _img_buf.getvalue()
            except Exception:
                _vaf_chart_bytes = None

            _vaf_mutant_cache = st.session_state.get("_vaf_mutant_assay_cache", "—")
            _vaf_wt_cache = st.session_state.get("_vaf_wt_assay_cache", "—")
            _vaf_summary_rows = [
                ["Parameter", "Value"],
                ["Mutant assay", _vaf_mutant_cache],
                ["Wild-type assay", _vaf_wt_cache],
                ["Samples analyzed", str(len(_vaf_results))],
                ["Samples with mutant detected", f"{_n_detected} / {len(_vaf_results)}"],
            ]
            _vaf_pdf_table_header = [_t['vaf_col_sample'], _t['vaf_col_lambda_mut'], _t['vaf_col_lambda_wt'],
                                      _t['vaf_col_fa'], _t['vaf_col_ci'], _t['vaf_col_conc_mut'], _t['vaf_col_detected']]
            _vaf_pdf_table_rows = [[str(v) for v in row.values()] for row in _vaf_rows]

            vaf_pdf_buffer = create_simple_pdf(
                report_title="AbsoluteGene",
                subtitle=_t['vaf_pdf_report'],
                description=_t['vaf_pdf_description'],
                summary_rows=_vaf_summary_rows,
                table_header=_vaf_pdf_table_header,
                table_rows=_vaf_pdf_table_rows,
                chart_png_bytes=_vaf_chart_bytes,
                chart_caption="Figure. Fractional Abundance (FA%, 95% CI) per sample. Red = mutant detected.",
                footer_note="AbsoluteGene — For research and educational use only. Not validated for clinical diagnostic purposes.",
                references=[
                    "Delta-method CI for ddPCR fractional abundance: Hindson CM et al. (2013). Nat Methods, 10(10), 1003-1005.",
                    "digital MIQE guidelines: Huggett JF et al. (2013). Clin Chem, 59(6), 892-902.",
                ]
            )
            st.download_button(
                "⬇️ " + _t['vaf_pdf_report'], data=vaf_pdf_buffer,
                file_name="vaf_report.pdf", mime="application/pdf", key="vaf_pdf_dl"
            )
    else:
        st.info(_t['vaf_no_data'])


# ═══════════════════════════════════════════════════════════════════════════════
# CLINICAL TOOLS TAB
# ═══════════════════════════════════════════════════════════════════════════════
with tab_clinical:
    st.markdown(f"### {_t['clinical_title']}")
    st.caption(_t['clinical_description'])
    st.markdown("---")

    clinical_mode = st.radio(
        _t['clinical_mode_label'],
        options=[_t['clinical_mode_mu'], _t['clinical_mode_rcv'],
                 _t['clinical_mode_precision'], _t['clinical_mode_comparison']],
        key="clinical_mode", horizontal=True
    )
    st.markdown("---")

    # ── MU Budget ──────────────────────────────────────────────────────────────
    if clinical_mode == _t['clinical_mode_mu']:
        st.caption(_t['mu_description'])

        _mu_gene_options = [f"{r['__gene__']} / {r['__group__']}" for r in data] if data else []
        mu_c1, mu_c2 = st.columns(2)
        with mu_c1:
            _mu_source_mode = st.radio(_t['mu_poisson_source'],
                                        options=[_t['mu_poisson_auto'], _t['mu_poisson_manual']],
                                        key="mu_source_mode")
        with mu_c2:
            if _mu_source_mode == _t['mu_poisson_auto'] and _mu_gene_options:
                _mu_selected = st.selectbox(_t['mu_gene_select'], options=_mu_gene_options, key="mu_gene_select")
                _mu_match = data[_mu_gene_options.index(_mu_selected)]
                _u_poisson_val = _mu_match.get("__conc_smp_cv__", 0.0)
                _u_poisson_val = 0.0 if (_u_poisson_val is None or np.isnan(_u_poisson_val)) else _u_poisson_val
                st.metric(_t['mu_poisson_source'], f"{_u_poisson_val:.2f}%")
            else:
                _u_poisson_val = st.number_input(_t['mu_poisson_source'], min_value=0.0, max_value=100.0,
                                                  value=5.0, step=0.1, key="mu_poisson_manual_val",
                                                  help=_t['mu_poisson_help'])

        mu_c3, mu_c4, mu_c5 = st.columns(3)
        with mu_c3:
            _u_pipetting_val = st.number_input(_t['mu_pipetting_label'], min_value=0.0, max_value=50.0,
                                                value=1.5, step=0.1, key="mu_pipetting_val",
                                                help=_t['mu_pipetting_help'])
        with mu_c4:
            _u_precision_val = st.number_input(_t['mu_precision_label'], min_value=0.0, max_value=50.0,
                                                value=0.0, step=0.1, key="mu_precision_val",
                                                help=_t['mu_precision_help'])
        with mu_c5:
            _mu_k = st.number_input(_t['mu_k_label'], min_value=1.0, max_value=3.0, value=2.0,
                                     step=0.5, key="mu_k_val")

        if st.button(_t['mu_calc_btn'], key="mu_calc_btn"):
            _mu_result = compute_mu_budget(_u_poisson_val, _u_pipetting_val, _u_precision_val, _mu_k)
            _val_label = _mu_selected if (_mu_source_mode == _t['mu_poisson_auto'] and _mu_gene_options) else "—"
            st.success(_t['mu_result_title'].format(value=_val_label, U=_mu_result["U_pct"], k=_mu_k))

            st.markdown(_t['mu_budget_table_title'])
            _mu_table_rows = [{_t['mu_col_source']: k, _t['mu_col_contribution']: f"{v:.3f}%"}
                               for k, v in _mu_result["components"].items()]
            _mu_table_rows.append({_t['mu_col_source']: _t['mu_col_combined'], _t['mu_col_contribution']: f"{_mu_result['u_c_pct']:.3f}%"})
            _mu_table_rows.append({_t['mu_col_source']: f"{_t['mu_col_expanded']} (k={_mu_k})", _t['mu_col_contribution']: f"{_mu_result['U_pct']:.3f}%"})
            st.dataframe(pd.DataFrame(_mu_table_rows), use_container_width=True)

            fig_mu = go.Figure(go.Bar(
                x=list(_mu_result["components"].keys()), y=list(_mu_result["components"].values()),
                marker_color="#00796b"
            ))
            fig_mu.update_layout(title="Uncertainty Component Contributions", yaxis_title="Relative uncertainty (%)", height=350)
            st.plotly_chart(fig_mu, use_container_width=True, key="mu_chart")

    # ── RCV ────────────────────────────────────────────────────────────────────
    elif clinical_mode == _t['clinical_mode_rcv']:
        st.caption(_t['rcv_description'])
        rcv_c1, rcv_c2 = st.columns(2)
        with rcv_c1:
            _rcv_result1 = st.number_input(_t['rcv_result1_label'], value=10.0, step=0.1, key="rcv_result1")
        with rcv_c2:
            _rcv_result2 = st.number_input(_t['rcv_result2_label'], value=12.0, step=0.1, key="rcv_result2")
        rcv_c3, rcv_c4, rcv_c5 = st.columns(3)
        with rcv_c3:
            _rcv_cv_a = st.number_input(_t['rcv_cv_analytical_label'], min_value=0.0, max_value=100.0,
                                         value=5.0, step=0.1, key="rcv_cv_a", help=_t['rcv_cv_analytical_help'])
        with rcv_c4:
            _rcv_cv_b = st.number_input(_t['rcv_cv_biological_label'], min_value=0.0, max_value=200.0,
                                         value=15.0, step=0.5, key="rcv_cv_b", help=_t['rcv_cv_biological_help'])
        with rcv_c5:
            _rcv_z = st.number_input(_t['rcv_z_label'], min_value=1.0, max_value=3.0, value=1.96,
                                      step=0.01, key="rcv_z")

        if st.button(_t['rcv_calc_btn'], key="rcv_calc_btn"):
            _rcv_res = compute_rcv(_rcv_result1, _rcv_result2, _rcv_cv_a, _rcv_cv_b, _rcv_z)
            if _rcv_res is None:
                st.error("Invalid input.")
            else:
                if _rcv_res["significant"]:
                    st.error(_t['rcv_result_significant'].format(change=_rcv_res["percent_change"], rcv=_rcv_res["rcv_pct"]))
                else:
                    st.success(_t['rcv_result_not_significant'].format(change=_rcv_res["percent_change"], rcv=_rcv_res["rcv_pct"]))
                rc1, rc2 = st.columns(2)
                rc1.metric("RCV (%)", f"±{_rcv_res['rcv_pct']:.2f}%")
                rc2.metric("Observed Change (%)", f"{_rcv_res['percent_change']:+.2f}%")

    # ── Precision Study ───────────────────────────────────────────────────────
    elif clinical_mode == _t['clinical_mode_precision']:
        st.caption(_t['precision_description'])
        _prec_n_days = st.number_input(_t['precision_n_days'], min_value=2, max_value=20, value=3, step=1, key="prec_n_days")
        _prec_day_arrays = []
        _prec_cols = st.columns(min(int(_prec_n_days), 4))
        for _d in range(int(_prec_n_days)):
            with _prec_cols[_d % len(_prec_cols)]:
                _txt = st.text_area(_t['precision_day_label'].format(i=_d + 1), key=f"prec_day_{_d}", height=120)
                _prec_day_arrays.append(parse_input_data(_txt))

        if st.button(_t['precision_calc_btn'], key="prec_calc_btn"):
            _prec_result = compute_precision_study(_prec_day_arrays)
            if _prec_result is None:
                st.error("Insufficient data — need at least 2 days with ≥2 replicates each.")
            else:
                pc1, pc2, pc3, pc4 = st.columns(4)
                pc1.metric(_t['precision_grand_mean'], f"{_prec_result['grand_mean']:.4f}")
                pc2.metric(_t['precision_repeatability_cv'], f"{_prec_result['repeatability_cv']:.2f}%")
                pc3.metric(_t['precision_between_day_cv'], f"{_prec_result['between_day_cv']:.2f}%")
                pc4.metric(_t['precision_total_cv'], f"{_prec_result['total_cv']:.2f}%")

                fig_prec = go.Figure()
                for _d, arr in enumerate(_prec_day_arrays):
                    if len(arr) > 0:
                        fig_prec.add_trace(go.Box(y=arr, name=f"Day {_d+1}", boxpoints="all"))
                fig_prec.update_layout(title="Precision Study — Per-Day Distribution", height=380)
                st.plotly_chart(fig_prec, use_container_width=True, key="prec_chart")

    # ── Method Comparison ─────────────────────────────────────────────────────
    else:
        st.caption(_t['comparison_description'])
        comp_c1, comp_c2 = st.columns(2)
        with comp_c1:
            _comp_m1_txt = st.text_area(_t['comparison_method1_label'], key="comp_m1", height=150)
        with comp_c2:
            _comp_m2_txt = st.text_area(_t['comparison_method2_label'], key="comp_m2", height=150)
        _comp_variance_ratio = st.number_input(_t['comparison_variance_ratio_label'], min_value=0.01, max_value=100.0,
                                                value=1.0, step=0.1, key="comp_var_ratio")

        if st.button(_t['comparison_calc_btn'], key="comp_calc_btn"):
            _m1_arr = parse_input_data(_comp_m1_txt)
            _m2_arr = parse_input_data(_comp_m2_txt)
            _ba_result = compute_bland_altman(_m1_arr, _m2_arr)
            _deming_result = compute_deming_regression(_m1_arr, _m2_arr, _comp_variance_ratio)

            if _ba_result is None:
                st.error("Need at least 2 paired values.")
            else:
                bc1, bc2, bc3 = st.columns(3)
                bc1.metric(_t['comparison_bias_label'], f"{_ba_result['mean_diff']:.4f}")
                bc2.metric(_t['comparison_loa_label'], f"{_ba_result['loa_low']:.3f} to {_ba_result['loa_high']:.3f}")
                if _deming_result:
                    bc3.metric(_t['comparison_deming_label'], f"y={_deming_result['slope']:.3f}x+{_deming_result['intercept']:.3f}")

                fig_ba = go.Figure()
                fig_ba.add_trace(go.Scatter(x=_ba_result["means"], y=_ba_result["diffs"], mode="markers",
                                             marker=dict(color="#00796b")))
                fig_ba.add_hline(y=_ba_result["mean_diff"], line_color="black", line_dash="solid")
                fig_ba.add_hline(y=_ba_result["loa_low"], line_color="red", line_dash="dash")
                fig_ba.add_hline(y=_ba_result["loa_high"], line_color="red", line_dash="dash")
                fig_ba.update_layout(title=_t['comparison_ba_chart_title'], xaxis_title="Mean of methods",
                                      yaxis_title="Difference (M2-M1)", height=380)
                st.plotly_chart(fig_ba, use_container_width=True, key="ba_chart")

                if _deming_result:
                    fig_dem = go.Figure()
                    fig_dem.add_trace(go.Scatter(x=_m1_arr, y=_m2_arr, mode="markers", marker=dict(color="#00796b"), name="Data"))
                    _x_line = np.array([min(_m1_arr), max(_m1_arr)])
                    _y_line = _deming_result["slope"] * _x_line + _deming_result["intercept"]
                    fig_dem.add_trace(go.Scatter(x=_x_line, y=_y_line, mode="lines", line=dict(color="red"), name="Deming fit"))
                    fig_dem.update_layout(title=_t['comparison_deming_chart_title'], xaxis_title="Method 1",
                                           yaxis_title="Method 2", height=380)
                    st.plotly_chart(fig_dem, use_container_width=True, key="deming_chart")


# ═══════════════════════════════════════════════════════════════════════════════
# CRM PRODUCTION TAB
# ═══════════════════════════════════════════════════════════════════════════════
with tab_crm:
    st.markdown(f"### {_t['crm_title']}")
    st.caption(_t['crm_description'])
    st.markdown("---")

    crm_mode = st.radio(
        _t['crm_mode_label'],
        options=[_t['crm_mode_homogeneity'], _t['crm_mode_stability'],
                 _t['crm_mode_uncertainty'], _t['crm_mode_coa']],
        key="crm_mode", horizontal=True
    )
    st.markdown("---")

    # ── Homogeneity Testing ───────────────────────────────────────────────────
    if crm_mode == _t['crm_mode_homogeneity']:
        st.caption(_t['homog_description'])
        _homog_n_units = st.number_input(_t['homog_n_units'], min_value=2, max_value=30, value=10, step=1, key="homog_n_units")
        _homog_unit_arrays = []
        _homog_cols = st.columns(min(int(_homog_n_units), 5))
        for _u in range(int(_homog_n_units)):
            with _homog_cols[_u % len(_homog_cols)]:
                _txt = st.text_area(_t['homog_unit_label'].format(i=_u + 1), key=f"homog_unit_{_u}", height=100)
                _homog_unit_arrays.append(parse_input_data(_txt))

        if st.button(_t['homog_calc_btn'], key="homog_calc_btn"):
            _homog_result = compute_homogeneity(_homog_unit_arrays)
            if _homog_result is None:
                st.error("Insufficient data — need at least 2 units with ≥2 replicates each.")
            else:
                st.session_state["_homog_result_cache"] = _homog_result
                if _homog_result["is_homogeneous"]:
                    st.success(_t['homog_result_homogeneous'].format(
                        F=_homog_result["F"], fcrit=_homog_result["F_crit"], p=_homog_result["p_value"]))
                else:
                    st.error(_t['homog_result_inhomogeneous'].format(
                        F=_homog_result["F"], fcrit=_homog_result["F_crit"], p=_homog_result["p_value"]))

                hc1, hc2, hc3 = st.columns(3)
                hc1.metric(_t['homog_grand_mean_label'], f"{_homog_result['grand_mean']:.5f}")
                hc2.metric(_t['homog_ubb_label'], f"{_homog_result['u_bb']:.5f}")
                hc3.metric("F / F_crit", f"{_homog_result['F']:.3f} / {_homog_result['F_crit']:.3f}")

                fig_homog = go.Figure()
                for _u, arr in enumerate(_homog_unit_arrays):
                    if len(arr) > 0:
                        fig_homog.add_trace(go.Box(y=arr, name=f"Unit {_u+1}", boxpoints="all"))
                fig_homog.add_hline(y=_homog_result["grand_mean"], line_dash="dash", line_color="black")
                fig_homog.update_layout(title="Homogeneity — Per-Unit Distribution", height=400)
                st.plotly_chart(fig_homog, use_container_width=True, key="homog_chart")

    # ── Stability Testing ─────────────────────────────────────────────────────
    elif crm_mode == _t['crm_mode_stability']:
        st.caption(_t['stab_description'])
        stab_c1, stab_c2 = st.columns(2)
        with stab_c1:
            _stab_time_txt = st.text_area(_t['stab_time_label'], key="stab_time", height=150)
        with stab_c2:
            _stab_value_txt = st.text_area(_t['stab_value_label'], key="stab_value", height=150)
        _stab_duration = st.number_input(_t['stab_duration_label'], min_value=0.1, value=365.0, step=1.0, key="stab_duration")

        if st.button(_t['stab_calc_btn'], key="stab_calc_btn"):
            _stab_time_arr = parse_input_data(_stab_time_txt)
            _stab_value_arr = parse_input_data(_stab_value_txt)
            _stab_result = compute_stability(_stab_time_arr, _stab_value_arr, _stab_duration)
            if _stab_result is None:
                st.error("Need at least 3 paired time/value points.")
            else:
                st.session_state["_stab_result_cache"] = _stab_result
                if _stab_result["is_stable"]:
                    st.success(_t['stab_result_stable'].format(p=_stab_result["p_value"]))
                else:
                    st.warning(_t['stab_result_unstable'].format(p=_stab_result["p_value"]))

                sc1, sc2, sc3 = st.columns(3)
                sc1.metric(_t['stab_slope_label'], f"{_stab_result['slope']:.6f}")
                sc2.metric(_t['stab_ustab_label'], f"{_stab_result['u_stab']:.6f}")
                sc3.metric("R²", f"{_stab_result['r_value']**2:.4f}")

                fig_stab = go.Figure()
                fig_stab.add_trace(go.Scatter(x=_stab_time_arr, y=_stab_value_arr, mode="markers",
                                               marker=dict(color="#00796b", size=9), name="Data"))
                _x_line = np.array([min(_stab_time_arr), max(_stab_time_arr)])
                _y_line = _stab_result["slope"] * _x_line + _stab_result["intercept"]
                fig_stab.add_trace(go.Scatter(x=_x_line, y=_y_line, mode="lines", line=dict(color="red"), name="Trend"))
                fig_stab.update_layout(title=_t['stab_chart_title'], xaxis_title="Time", yaxis_title="Value", height=400)
                st.plotly_chart(fig_stab, use_container_width=True, key="stab_chart")

    # ── Assigned Value & Uncertainty Budget ───────────────────────────────────
    elif crm_mode == _t['crm_mode_uncertainty']:
        st.caption(_t['unc_description'])

        _cached_homog = st.session_state.get("_homog_result_cache")
        _cached_stab = st.session_state.get("_stab_result_cache")
        _use_cached = st.checkbox(_t['unc_use_cached'], value=True, key="unc_use_cached")
        if _use_cached and not (_cached_homog or _cached_stab):
            st.info(_t['unc_no_cache'])

        unc_c1, unc_c2 = st.columns(2)
        with unc_c1:
            _unc_assigned_value = st.number_input(_t['unc_assigned_value_label'], value=100.0, step=0.1, key="unc_assigned_value")
            _unc_u_char = st.number_input(_t['unc_u_char_label'], min_value=0.0, value=1.0, step=0.1,
                                           key="unc_u_char", help=_t['unc_u_char_help'])
        with unc_c2:
            _default_u_bb = _cached_homog["u_bb"] if (_use_cached and _cached_homog) else 0.0
            _unc_u_bb = st.number_input(_t['unc_u_bb_label'], min_value=0.0, value=float(_default_u_bb), step=0.01, key="unc_u_bb")
            _default_u_stab = _cached_stab["u_stab"] if (_use_cached and _cached_stab) else 0.0
            _unc_u_stab = st.number_input(_t['unc_u_stab_label'], min_value=0.0, value=float(_default_u_stab), step=0.01, key="unc_u_stab")

        _unc_k = st.number_input(_t['unc_k_label'], min_value=1.0, max_value=3.0, value=2.0, step=0.5, key="unc_k")

        if st.button(_t['unc_calc_btn'], key="unc_calc_btn"):
            _unc_result = compute_assigned_value_uncertainty(_unc_assigned_value, _unc_u_char, _unc_u_bb, _unc_u_stab, k=_unc_k)
            st.session_state["_unc_result_cache"] = _unc_result
            st.success(_t['unc_result_title'].format(
                value=_unc_result["assigned_value"], U=_unc_result["U"], k=_unc_k, urel=_unc_result["U_rel_pct"]))

            _unc_table_rows = [{_t['mu_col_source']: k, _t['mu_col_contribution']: f"{v:.5f}"}
                                for k, v in _unc_result["components"].items()]
            _unc_table_rows.append({_t['mu_col_source']: _t['mu_col_combined'], _t['mu_col_contribution']: f"{_unc_result['u_c']:.5f}"})
            _unc_table_rows.append({_t['mu_col_source']: f"{_t['mu_col_expanded']} (k={_unc_k})", _t['mu_col_contribution']: f"{_unc_result['U']:.5f}"})
            st.dataframe(pd.DataFrame(_unc_table_rows), use_container_width=True)

            fig_unc = go.Figure(go.Bar(
                x=list(_unc_result["components"].keys()), y=list(_unc_result["components"].values()),
                marker_color="#00695c"
            ))
            fig_unc.update_layout(title="Uncertainty Component Contributions (absolute units)", height=350)
            st.plotly_chart(fig_unc, use_container_width=True, key="unc_chart")

    # ── CoA Generator ──────────────────────────────────────────────────────────
    else:
        st.caption(_t['coa_description'])
        _cached_unc = st.session_state.get("_unc_result_cache")

        coa_c1, coa_c2 = st.columns(2)
        with coa_c1:
            _coa_material = st.text_input(_t['coa_material_name'], value="", key="coa_material")
            _coa_lot = st.text_input(_t['coa_lot_number'], value="", key="coa_lot")
            _coa_producer = st.text_input(_t['coa_producer'], value="", key="coa_producer")
        with coa_c2:
            _default_val = _cached_unc["assigned_value"] if _cached_unc else 100.0
            _default_U = _cached_unc["U"] if _cached_unc else 2.0
            _default_k = _cached_unc["k"] if _cached_unc else 2.0
            _coa_value = st.number_input(_t['coa_assigned_value_label'], value=float(_default_val), step=0.01, key="coa_value")
            _coa_unit = st.text_input(_t['coa_unit_label'], value="copies/µL", key="coa_unit")
            _coa_U = st.number_input(_t['coa_expanded_unc_label'], value=float(_default_U), step=0.01, key="coa_U")
            _coa_k = st.number_input(_t['coa_k_label'], value=float(_default_k), step=0.5, key="coa_k")

        _coa_traceability = st.text_area(_t['coa_traceability_label'], value=_t['coa_traceability_default'], key="coa_traceability")
        _coa_validity = st.text_input(_t['coa_validity_label'], value="", key="coa_validity")

        if st.button(_t['coa_generate_btn'], key="coa_generate_btn"):
            _coa_summary_rows = [
                ["Parameter", "Value"],
                ["Material", _coa_material or "—"],
                ["Batch/Lot", _coa_lot or "—"],
                ["Producer/Laboratory", _coa_producer or "—"],
                ["Assigned Value", f"{_coa_value} {_coa_unit}"],
                ["Expanded Uncertainty (U)", f"± {_coa_U} {_coa_unit} (k={_coa_k})"],
                ["Validity / Shelf Life", _coa_validity or "—"],
            ]
            coa_pdf_buffer = create_simple_pdf(
                report_title="Certificate of Analysis",
                subtitle=_coa_material or "AbsoluteGene CRM Certificate",
                description=_coa_traceability,
                summary_rows=_coa_summary_rows,
                table_header=None, table_rows=None, chart_png_bytes=None, chart_caption=None,
                footer_note=("This certificate was generated by AbsoluteGene for research and internal "
                              "documentation purposes. It does not constitute a formally accredited "
                              "Certificate of Analysis under ISO 17034 unless issued by an accredited "
                              "reference material producer following full validation."),
                references=[
                    "ISO Guide 35:2017. Reference materials — Guidance for characterization and assessment of homogeneity and stability.",
                    "ISO 17034:2016. General requirements for the competence of reference material producers.",
                    "Linsinger TP et al. (2001). Homogeneity and stability of reference materials. Accred Qual Assur, 6, 20-25.",
                ]
            )
            st.download_button(
                _t['coa_download_btn'], data=coa_pdf_buffer,
                file_name="certificate_of_analysis.pdf", mime="application/pdf", key="coa_pdf_dl"
            )


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
