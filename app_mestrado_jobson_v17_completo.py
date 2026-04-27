# -*- coding: utf-8 -*-
"""
Dashboard Executivo & Estatístico - Stewardship Antimicrobiano
Mestrado Fiocruz - Versão expandida V18

Esta versão mantém a planilha/formulário já existentes e cria novas análises
apenas a partir dos campos já disponíveis na aba PRESCRIÇÕES.

Como usar:
1) Coloque este arquivo .py na mesma pasta da planilha:
   Sistema de Monitoramento - Mestrado Jobson (V14.0).xlsx
2) Execute:
   python -m streamlit run app_mestrado_jobson_v18_expandido.py

Dependências principais:
pip install streamlit pandas numpy scipy statsmodels plotly openpyxl
"""

import itertools
import re
import unicodedata
from pathlib import Path

import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import scipy.stats as stats
import statsmodels.api as sm
import streamlit as st

# ============================================================
# 1. CONFIGURAÇÃO GERAL
# ============================================================

st.set_page_config(
    page_title="Stewardship Antimicrobiano - Mestrado Fiocruz",
    page_icon="💊",
    layout="wide",
)

px.defaults.template = "plotly_white"
px.defaults.color_continuous_scale = "Viridis"

st.title("📊 Dashboard Executivo & Estatístico - Stewardship Antimicrobiano")
st.caption(
    "V18 expandida: mantém a planilha/formulário existentes e adiciona Tabela 1, "
    "auditoria de qualidade, AWaRe aproximado, ASI, profilaxia cirúrgica, sobreposição potencial, "
    "ranking de antimicrobianos, bootstrap e simulação epidemiológica."
)
st.markdown("---")

ARQUIVO_PADRAO = "Sistema de Monitoramento - Mestrado Jobson (V14.0).xlsx"
ABA_PADRAO = "PRESCRIÇÕES"
ALFA = 0.05
RANDOM_SEED = 42

COLUNAS_NUMERICAS = [
    "CUSTO_TOTAL_R$", "Idade", "DOT_Exato", "Dose_diaria", "DDD_OMS",
    "DDD_TOTAL", "TEMPO_PROTOCOLO", "Tempo_Prescrito", "CUSTO_MG_R$"
]

HIERARQUIA_DESFECHO = {
    "Óbito": 5,
    "Obito": 5,
    "Transferido": 4,
    "Transferência": 4,
    "Transferencia": 4,
    "EVASAO": 3,
    "Evasão": 3,
    "Evasao": 3,
    "Em tratamento (Ativo)": 2,
    "Alta hospitalar": 1,
}

# Mapeamento aproximado AWaRe para os antimicrobianos mais prováveis na base.
# Recomenda-se revisar/validar antes de usar em artigo, conforme lista OMS vigente.
AWARE_MAP = {
    "amoxicilina": "Access",
    "amoxicilina clavulanato": "Access",
    "ampicilina": "Access",
    "ampicilina sulbactam": "Access",
    "cefazolina": "Access",
    "cefalexina": "Access",
    "clindamicina": "Access",
    "gentamicina": "Access",
    "metronidazol": "Access",
    "nitrofurantoina": "Access",
    "sulfametoxazol trimetoprim": "Access",
    "trimetoprim sulfametoxazol": "Access",
    "ceftriaxona": "Watch",
    "cefotaxima": "Watch",
    "cefepime": "Watch",
    "ceftazidima": "Watch",
    "ciprofloxacino": "Watch",
    "levofloxacino": "Watch",
    "moxifloxacino": "Watch",
    "piperacilina tazobactam": "Watch",
    "azitromicina": "Watch",
    "claritromicina": "Watch",
    "vancomicina": "Watch",
    "teicoplanina": "Watch",
    "meropenem": "Reserve",
    "imipenem": "Reserve",
    "ertapenem": "Reserve",
    "linezolida": "Reserve",
    "linezolid": "Reserve",
    "daptomicina": "Reserve",
    "polimixina": "Reserve",
    "polimixina b": "Reserve",
    "colistina": "Reserve",
    "tigeciclina": "Reserve",
    "ceftazidima avibactam": "Reserve",
}

# ============================================================
# 2. UTILIDADES TEXTUAIS E ESTATÍSTICAS
# ============================================================


def remover_acentos(texto):
    if pd.isna(texto):
        return ""
    texto = str(texto).strip()
    return "".join(
        ch for ch in unicodedata.normalize("NFKD", texto)
        if not unicodedata.combining(ch)
    )


def norm_text(texto):
    x = remover_acentos(texto).lower()
    x = re.sub(r"[^a-z0-9]+", " ", x)
    return re.sub(r"\s+", " ", x).strip()


def padronizar_setor(valor):
    x = norm_text(valor).upper()
    if "PRE" in x and ("OP" in x or "OPER" in x):
        return "Pré-operatório"
    if "OBST" in x:
        return "Enfermaria Obstétrica"
    if "CIR" in x:
        return "Enfermaria Cirúrgica"
    if "CLIN" in x:
        return "Enfermaria Clínica"
    if "UTI" in x:
        return "UTI"
    return str(valor).strip() if pd.notna(valor) else np.nan


def excel_serial_para_data(serie):
    return pd.to_datetime(serie, unit="D", origin="1899-12-30", errors="coerce")


def converter_data_mista(serie):
    """Converte datas vindas como serial do Excel ou texto."""
    numerica = pd.to_numeric(serie, errors="coerce")
    out = pd.Series(pd.NaT, index=serie.index, dtype="datetime64[ns]")
    mask_num = numerica.notna()
    if mask_num.any():
        out.loc[mask_num] = excel_serial_para_data(numerica.loc[mask_num])
    mask_txt = ~mask_num & serie.notna()
    if mask_txt.any():
        out.loc[mask_txt] = pd.to_datetime(serie.loc[mask_txt], errors="coerce", dayfirst=True)
    return out


def classificar_ddd(ratio):
    if pd.isna(ratio) or ratio <= 0:
        return "Não Avaliado"
    if 0.9 <= ratio <= 1.1:
        return "Adequada à DDD (±10%)"
    if ratio < 0.9:
        return "Subdose"
    return "Superdose"


def wilson_ic(x, n, alpha=0.05):
    if n == 0 or pd.isna(n):
        return np.nan, np.nan
    z = stats.norm.ppf(1 - alpha / 2)
    p = x / n
    den = 1 + z**2 / n
    centro = (p + z**2 / (2 * n)) / den
    margem = z * np.sqrt((p * (1 - p) + z**2 / (4 * n)) / n) / den
    return centro - margem, centro + margem


def formatar_pct_ic(x, n):
    if n == 0 or pd.isna(n):
        return "NA"
    li, ls = wilson_ic(x, n)
    return f"{100*x/n:.1f}% (IC95% {100*li:.1f}–{100*ls:.1f})"


def formatar_p(p):
    if p is None or pd.isna(p):
        return "NA"
    if p < 0.001:
        return "<0,001"
    return f"{p:.4f}".replace(".", ",")


def interpretar_p(p, alfa=ALFA):
    if p is None or pd.isna(p):
        return "Não interpretável por ausência ou insuficiência de dados."
    if p < alfa:
        return "Associação/diferença estatisticamente significativa ao nível de 5%."
    return "Sem evidência estatística suficiente de associação/diferença ao nível de 5%; isso não prova ausência de efeito."


def classe_decisao_p(p, alfa=ALFA):
    if p is None or pd.isna(p):
        return "Indeterminado"
    return "Significativo" if p < alfa else "Não significativo"


def interpretar_rho(rho):
    if rho is None or pd.isna(rho):
        return "Correlação não estimável."
    abs_r = abs(rho)
    direcao = "positiva" if rho > 0 else "negativa" if rho < 0 else "nula"
    if abs_r < 0.20:
        intensidade = "muito fraca"
    elif abs_r < 0.40:
        intensidade = "fraca"
    elif abs_r < 0.60:
        intensidade = "moderada"
    elif abs_r < 0.80:
        intensidade = "forte"
    else:
        intensidade = "muito forte"
    return f"Correlação {direcao} {intensidade}."


def mostrar_resultado_teste(nome, estatistica=None, p=None, detalhes="", alfa=ALFA):
    partes = [f"**{nome}**"]
    if estatistica is not None and pd.notna(estatistica):
        partes.append(f"estatística={estatistica:.3f}".replace(".", ","))
    if p is not None and pd.notna(p):
        partes.append(f"p-valor={formatar_p(p)}")
    if detalhes:
        partes.append(detalhes)
    st.info("; ".join(partes) + ".")
    interp = interpretar_p(p, alfa)
    if p is not None and pd.notna(p) and p < alfa:
        st.success(interp)
    else:
        st.warning(interp)


def descrever_num(serie):
    s = pd.to_numeric(serie, errors="coerce").dropna()
    if len(s) == 0:
        return pd.Series({
            "n": 0, "média": np.nan, "DP": np.nan, "mediana": np.nan,
            "Q1": np.nan, "Q3": np.nan, "mín": np.nan, "máx": np.nan,
        })
    return pd.Series({
        "n": int(len(s)),
        "média": s.mean(),
        "DP": s.std(ddof=1),
        "mediana": s.median(),
        "Q1": s.quantile(0.25),
        "Q3": s.quantile(0.75),
        "mín": s.min(),
        "máx": s.max(),
    })


def desfecho_final(grupo):
    vals = grupo.dropna().astype(str).unique().tolist()
    if not vals:
        return np.nan
    return sorted(vals, key=lambda v: HIERARQUIA_DESFECHO.get(v, 0), reverse=True)[0]


def moda_segura(serie):
    m = serie.dropna().mode()
    return m.iloc[0] if len(m) else np.nan


def nome_coluna_prescricao(df):
    for c in ["ID_PRESCRIÇÃO", "ID_PRESCRICAO", "ID_Prescricao", "ID_Prescrição"]:
        if c in df.columns:
            return c
    return None


def classificar_aware(atb):
    txt = norm_text(atb)
    if not txt:
        return "Não classificado"
    # prioriza chaves longas para evitar classificar piperacilina/tazo só por termo parcial
    for chave in sorted(AWARE_MAP, key=len, reverse=True):
        if chave in txt:
            return AWARE_MAP[chave]
    return "Não classificado"


def calcular_asi(row):
    """Antimicrobial Spectrum Index simples, derivado de ESPECTRO/AWaRe."""
    espectro = norm_text(row.get("ESPECTRO", ""))
    aware = row.get("AWaRe", "Não classificado")
    if aware == "Reserve":
        return 4
    if "reserva" in espectro or "critico" in espectro:
        return 4
    if "amplo" in espectro:
        return 3
    if "moder" in espectro or "intermedi" in espectro:
        return 2
    if "estreito" in espectro or "reduz" in espectro:
        return 1
    if aware == "Watch":
        return 3
    if aware == "Access":
        return 1
    return np.nan


def detectar_profilaxia(row):
    campos = " ".join(str(row.get(c, "")) for c in ["Motivo", "Sindrome", "AUD_PRESCRICAO"])
    txt = norm_text(campos)
    return int(("profil" in txt) or ("cirurg" in txt and "infecc" not in txt))


def otimizar_memoria(df):
    df = df.copy()
    for col in df.select_dtypes(include=["object"]).columns:
        nunique = df[col].nunique(dropna=True)
        if 0 < nunique <= max(50, len(df) * 0.5):
            df[col] = df[col].astype("category")
    return df

# ============================================================
# 3. CARREGAMENTO, PREPARAÇÃO E ENGENHARIA DE VARIÁVEIS
# ============================================================


@st.cache_data(show_spinner=False)
def carregar_arquivo_cache(caminho_ou_bytes, uploaded_name=None):
    """Carrega Excel/CSV em cache. Aceita caminho local ou bytes de upload."""
    if uploaded_name:
        nome = uploaded_name.lower()
        if nome.endswith((".xlsx", ".xls")):
            xls = pd.ExcelFile(caminho_ou_bytes)
            aba = ABA_PADRAO if ABA_PADRAO in xls.sheet_names else xls.sheet_names[0]
            return pd.read_excel(caminho_ou_bytes, sheet_name=aba)
        return pd.read_csv(caminho_ou_bytes, sep=";", decimal=",")

    caminho = Path(caminho_ou_bytes)
    if caminho.suffix.lower() in [".xlsx", ".xls"]:
        xls = pd.ExcelFile(caminho)
        aba = ABA_PADRAO if ABA_PADRAO in xls.sheet_names else xls.sheet_names[0]
        return pd.read_excel(caminho, sheet_name=aba)
    return pd.read_csv(caminho, sep=";", decimal=",")


def carregar_arquivo(uploaded_file=None):
    if uploaded_file is not None:
        return carregar_arquivo_cache(uploaded_file, uploaded_file.name)
    caminho_excel = Path(ARQUIVO_PADRAO)
    if caminho_excel.exists():
        return carregar_arquivo_cache(str(caminho_excel))
    caminho_csv = Path("PRESCRIÇÕES.csv")
    if caminho_csv.exists():
        return carregar_arquivo_cache(str(caminho_csv))
    raise FileNotFoundError(
        f"Coloque o arquivo '{ARQUIVO_PADRAO}' na mesma pasta do app, "
        "ou envie a planilha pela barra lateral."
    )


def preparar_variaveis(df):
    df = df.copy()
    df.columns = [str(c).strip() for c in df.columns]
    if "ID_Internacao" not in df.columns:
        raise ValueError("A coluna 'ID_Internacao' não foi encontrada.")

    df = df.dropna(subset=["ID_Internacao"]).copy()

    for col in COLUNAS_NUMERICAS:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")

    for col in ["Carimbo", "Inicio", "Fim", "Data", "Data_Prescricao"]:
        if col in df.columns:
            df[col + "_Data"] = converter_data_mista(df[col])

    df["Setor_Padronizado"] = df["Setor"].apply(padronizar_setor) if "Setor" in df.columns else np.nan

    if "Idade" in df.columns:
        bins = [0, 17, 39, 59, 130]
        labels = ["1. Pediátrico (<18)", "2. Adulto jovem (18–39)", "3. Adulto (40–59)", "4. Idoso (60+)"]
        df["Faixa_Etaria"] = pd.cut(df["Idade"], bins=bins, labels=labels, right=True)

    if {"Dose_diaria", "DDD_OMS"}.issubset(df.columns):
        df["Ratio_DDD"] = np.where(df["DDD_OMS"] > 0, (df["Dose_diaria"] / 1000) / df["DDD_OMS"], np.nan)
        df["Status_DDD"] = df["Ratio_DDD"].apply(classificar_ddd)
        df["DDD_Adequada"] = df["Status_DDD"].eq("Adequada à DDD (±10%)").astype(int)
        df["DDD_Inadequada"] = df["Status_DDD"].isin(["Subdose", "Superdose"]).astype(int)
        df["Alerta_Outlier_Dose"] = np.where((df["Ratio_DDD"] >= 1.5) | ((df["Ratio_DDD"] > 0) & (df["Ratio_DDD"] <= 0.5)), 1, 0)
    else:
        df["Ratio_DDD"] = np.nan
        df["Status_DDD"] = "Não Avaliado"
        df["DDD_Adequada"] = 0
        df["DDD_Inadequada"] = 0
        df["Alerta_Outlier_Dose"] = 0

    # Colunas garantidas para o restante do app
    for col in ["CUSTO_TOTAL_R$", "DOT_Exato", "Idade", "DDD_TOTAL", "Tempo_Prescrito"]:
        if col not in df.columns:
            df[col] = np.nan
    for col in ["ATB", "CLASSE", "Sindrome", "Motivo", "Sexo", "Desfecho", "AUD_PRESCRICAO", "AUD_ADMINISTRACAO", "Cultura", "ESPECTRO", "Determinou_Tempo?"]:
        if col not in df.columns:
            df[col] = np.nan

    aud_presc = df["AUD_PRESCRICAO"].astype(str)
    aud_admin = df["AUD_ADMINISTRACAO"].astype(str)
    desfecho = df["Desfecho"].astype(str)
    cultura = df["Cultura"].astype(str)
    espectro = df["ESPECTRO"].astype(str)
    tempo_det = df["Determinou_Tempo?"].astype(str)

    df["Prescricao_Nao_Informada"] = aud_presc.str.contains("não informado|nao informado", case=False, na=False).astype(int)
    df["Prescricao_Avaliavel"] = (~aud_presc.str.contains("não informado|nao informado|individualizada", case=False, na=False)).astype(int)
    df["Erro_Prescricao_Avaliavel"] = np.where(
        df["Prescricao_Avaliavel"].eq(1),
        aud_presc.str.contains("prolongada|curta|acima|abaixo", case=False, na=False).astype(int),
        np.nan,
    )

    df["Erro_Administracao"] = aud_admin.str.contains("incompleta|suspensão precoce|suspensao precoce", case=False, na=False).astype(int)
    df["Administracao_Prolongada"] = aud_admin.str.contains("prolongada", case=False, na=False).astype(int)
    df["Administracao_Completa"] = aud_admin.str.contains("ciclo completo", case=False, na=False).astype(int)
    df["Administracao_Sem_Base"] = aud_admin.str.contains("sem base", case=False, na=False).astype(int)

    df["Obito"] = desfecho.str.contains("óbito|obito", case=False, na=False).astype(int)
    df["Cultura_Sim"] = cultura.str.contains("^sim|sim", case=False, na=False).astype(int)
    df["Amplo_Espectro"] = espectro.str.contains("amplo", case=False, na=False).astype(int)
    df["Tempo_Determinado"] = tempo_det.str.contains("^sim", case=False, na=False).astype(int)

    # Novas variáveis derivadas sem alterar a planilha
    df["AWaRe"] = df["ATB"].apply(classificar_aware)
    df["ASI"] = df.apply(calcular_asi, axis=1)
    df["Profilaxia_Cirurgica"] = df.apply(detectar_profilaxia, axis=1)
    tempo_ref = pd.to_numeric(df["Tempo_Prescrito"], errors="coerce")
    tempo_ref = tempo_ref.fillna(pd.to_numeric(df["DOT_Exato"], errors="coerce"))
    df["Tempo_Analise_Profilaxia"] = tempo_ref
    df["Profilaxia_Maior_24h"] = np.where((df["Profilaxia_Cirurgica"] == 1) & (tempo_ref > 1), 1, 0)
    df["Profilaxia_Maior_48h"] = np.where((df["Profilaxia_Cirurgica"] == 1) & (tempo_ref > 2), 1, 0)
    df["Custo_Inconformidade"] = np.where(
        (df["Erro_Administracao"].eq(1)) | (df["DDD_Inadequada"].eq(1)) | (df["Profilaxia_Maior_24h"].eq(1)) | (df["Alerta_Outlier_Dose"].eq(1)),
        pd.to_numeric(df["CUSTO_TOTAL_R$"], errors="coerce").fillna(0),
        0,
    )
    return otimizar_memoria(df)


def detectar_sobreposicao_potencial(df):
    """Detecta possível double coverage sem alterar a planilha.

    Se houver Inicio_Data e Fim_Data, busca intervalos sobrepostos por internação.
    Se não houver datas, sinaliza internações com >=2 ATB diferentes e mesmo ESPECTRO/CLASSE.
    """
    if df.empty:
        return pd.DataFrame()
    linhas = []
    tem_datas = "Inicio_Data" in df.columns and "Fim_Data" in df.columns and df["Inicio_Data"].notna().any() and df["Fim_Data"].notna().any()
    for id_int, g in df.groupby("ID_Internacao", dropna=True):
        gg = g.copy()
        gg["_atb"] = gg["ATB"].astype(str)
        if gg["_atb"].nunique() < 2:
            continue
        for _, a in gg.iterrows():
            for _, b in gg.iterrows():
                if a.name >= b.name or a["_atb"] == b["_atb"]:
                    continue
                mesma_classe = str(a.get("CLASSE", "")) == str(b.get("CLASSE", "")) and str(a.get("CLASSE", "")) not in ["", "nan", "NaN"]
                mesmo_espectro = str(a.get("ESPECTRO", "")) == str(b.get("ESPECTRO", "")) and str(a.get("ESPECTRO", "")) not in ["", "nan", "NaN"]
                if not (mesma_classe or mesmo_espectro):
                    continue
                sobrepoe = True
                metodo = "potencial sem datas"
                if tem_datas and pd.notna(a.get("Inicio_Data")) and pd.notna(a.get("Fim_Data")) and pd.notna(b.get("Inicio_Data")) and pd.notna(b.get("Fim_Data")):
                    sobrepoe = (a["Inicio_Data"] <= b["Fim_Data"]) and (b["Inicio_Data"] <= a["Fim_Data"])
                    metodo = "intervalo temporal" if sobrepoe else "sem sobreposição temporal"
                if sobrepoe:
                    linhas.append({
                        "ID_Internacao": id_int,
                        "ATB 1": a["_atb"],
                        "ATB 2": b["_atb"],
                        "Classe semelhante": mesma_classe,
                        "Espectro semelhante": mesmo_espectro,
                        "Método": metodo,
                    })
    return pd.DataFrame(linhas).drop_duplicates() if linhas else pd.DataFrame(columns=["ID_Internacao", "ATB 1", "ATB 2", "Classe semelhante", "Espectro semelhante", "Método"])


def agregar_internacoes(df):
    col_presc = nome_coluna_prescricao(df)
    presc_agg = (col_presc, "nunique") if col_presc else ("ID_Internacao", "size")
    agg = df.groupby("ID_Internacao", dropna=True).agg(
        n_prescricoes=presc_agg,
        idade=("Idade", "first"),
        sexo=("Sexo", "first"),
        setor_principal=("Setor_Padronizado", moda_segura),
        custo_total=("CUSTO_TOTAL_R$", "sum"),
        custo_inconformidade=("Custo_Inconformidade", "sum"),
        dot_total=("DOT_Exato", "sum"),
        ddd_total=("DDD_TOTAL", "sum"),
        n_atb=("ATB", "nunique"),
        amplo_espectro=("Amplo_Espectro", "max"),
        asi_total=("ASI", "sum"),
        asi_medio=("ASI", "mean"),
        aware_watch_reserve=("AWaRe", lambda x: int(any(v in ["Watch", "Reserve"] for v in x.astype(str)))),
        erro_admin_algum=("Erro_Administracao", "max"),
        tempo_determinado_algum=("Tempo_Determinado", "max"),
        cultura_alguma=("Cultura_Sim", "max"),
        obito=("Obito", "max"),
        prof_cirurgica_alguma=("Profilaxia_Cirurgica", "max"),
        prof_maior_24h_alguma=("Profilaxia_Maior_24h", "max"),
        outlier_dose_algum=("Alerta_Outlier_Dose", "max"),
        desfecho_final=("Desfecho", desfecho_final),
        desfechos_registrados=("Desfecho", lambda x: " | ".join(pd.unique(x.dropna().astype(str)))),
    ).reset_index()
    agg["multiplos_desfechos"] = agg["desfechos_registrados"].astype(str).str.contains(r"\|").astype(int)
    return otimizar_memoria(agg)

# ============================================================
# 4. TESTES, TABELAS E BOOTSTRAP
# ============================================================


def teste_categorico(tabela):
    tabela = tabela.copy()
    if tabela.empty or tabela.to_numpy().sum() == 0:
        return "Não aplicável", np.nan, np.nan, np.nan, np.nan
    if tabela.shape == (2, 2):
        odds, p = stats.fisher_exact(tabela.to_numpy())
        return "Fisher exato", p, odds, np.nan, np.nan
    chi2, p, dof, exp = stats.chi2_contingency(tabela)
    n = tabela.to_numpy().sum()
    k = min(tabela.shape) - 1
    v_cramer = np.sqrt(chi2 / (n * k)) if n > 0 and k > 0 else np.nan
    menor_esperado = np.min(exp) if exp.size else np.nan
    return "Qui-quadrado", p, v_cramer, menor_esperado, chi2


def tabela_2x2_binaria(df, exposicao, desfecho):
    tab = pd.crosstab(df[exposicao], df[desfecho])
    return tab.reindex(index=[0, 1], columns=[0, 1], fill_value=0)


def fisher_or_com_ic(tabela):
    arr = tabela.to_numpy()
    odds, p = stats.fisher_exact(arr)
    try:
        t22 = sm.stats.Table2x2(arr, shift_zeros=True)
        li, ls = t22.oddsratio_confint()
    except Exception:
        li, ls = np.nan, np.nan
    return odds, li, ls, p


def kruskal_por_grupo(df, grupo, valor):
    dados = []
    for nome, g in df.groupby(grupo, dropna=True):
        s = pd.to_numeric(g[valor], errors="coerce").dropna()
        if len(s) > 0:
            dados.append((nome, s))
    if len(dados) < 2:
        return None
    stat, p = stats.kruskal(*[s for _, s in dados])
    return stat, p, dados


def pares_mannwhitney_holm(dados):
    linhas = []
    for (g1, s1), (g2, s2) in itertools.combinations(dados, 2):
        if len(s1) > 0 and len(s2) > 0:
            u, p = stats.mannwhitneyu(s1, s2, alternative="two-sided")
            linhas.append({"Grupo 1": g1, "Grupo 2": g2, "U": u, "p bruto": p})
    if not linhas:
        return pd.DataFrame()
    out = pd.DataFrame(linhas).sort_values("p bruto").reset_index(drop=True)
    m = len(out)
    holm = []
    maior_ajustado = 0
    for i, p_val in enumerate(out["p bruto"]):
        ajustado = min((m - i) * p_val, 1)
        maior_ajustado = max(maior_ajustado, ajustado)
        holm.append(maior_ajustado)
    out["p Holm"] = holm
    out["Decisão"] = out["p Holm"].apply(classe_decisao_p)
    out["Interpretação"] = out["p Holm"].apply(interpretar_p)
    return out


def testar_kruskal_e_posthoc(df, grupo, valor, titulo):
    res = kruskal_por_grupo(df, grupo, valor)
    if res is None:
        st.warning(f"Dados insuficientes para {titulo}.")
        return pd.DataFrame()
    h, p, dados = res
    mostrar_resultado_teste(titulo, estatistica=h, p=p)
    posthoc = pd.DataFrame()
    if p < ALFA:
        st.write("**Pós-teste:** Mann-Whitney par-a-par com correção de Holm-Bonferroni.")
        posthoc = pares_mannwhitney_holm(dados)
        st.dataframe(posthoc, use_container_width=True)
    return posthoc


def correlacoes_spearman(df, variaveis):
    linhas = []
    for v1, v2 in itertools.combinations(variaveis, 2):
        base = df[[v1, v2]].dropna()
        if len(base) >= 5 and base[v1].nunique() > 1 and base[v2].nunique() > 1:
            rho, p = stats.spearmanr(base[v1], base[v2])
        else:
            rho, p = np.nan, np.nan
        linhas.append({
            "Variável 1": v1,
            "Variável 2": v2,
            "rho Spearman": rho,
            "p-valor": p,
            "p formatado": formatar_p(p),
            "Decisão": classe_decisao_p(p),
            "Interpretação estatística": interpretar_p(p),
            "Interpretação da força": interpretar_rho(rho),
        })
    return pd.DataFrame(linhas)


def bootstrap_ci_median(serie, n_boot=2000, alpha=0.05, seed=RANDOM_SEED):
    x = pd.to_numeric(serie, errors="coerce").dropna().to_numpy()
    if len(x) < 2:
        return np.nan, np.nan, np.nan
    rng = np.random.default_rng(seed)
    meds = np.array([np.median(rng.choice(x, size=len(x), replace=True)) for _ in range(n_boot)])
    return np.median(x), np.quantile(meds, alpha/2), np.quantile(meds, 1-alpha/2)


def bootstrap_ci_spearman(df, xcol, ycol, n_boot=2000, alpha=0.05, seed=RANDOM_SEED):
    base = df[[xcol, ycol]].dropna()
    if len(base) < 5 or base[xcol].nunique() < 2 or base[ycol].nunique() < 2:
        return np.nan, np.nan, np.nan
    rho_obs, _ = stats.spearmanr(base[xcol], base[ycol])
    rng = np.random.default_rng(seed)
    vals = []
    arr = base.to_numpy()
    n = len(arr)
    for _ in range(n_boot):
        sample = arr[rng.integers(0, n, size=n)]
        if len(np.unique(sample[:, 0])) > 1 and len(np.unique(sample[:, 1])) > 1:
            vals.append(stats.spearmanr(sample[:, 0], sample[:, 1])[0])
    if not vals:
        return rho_obs, np.nan, np.nan
    return rho_obs, np.quantile(vals, alpha/2), np.quantile(vals, 1-alpha/2)


def gerar_tabela1(df, outcome_col, variaveis):
    if outcome_col not in df.columns or df[outcome_col].nunique(dropna=True) < 2:
        return pd.DataFrame()
    vc = df[outcome_col].value_counts(dropna=False)
    n0 = int(vc.get(0, 0))
    n1 = int(vc.get(1, 0))
    col0 = f"Ausente/Não (n={n0})"
    col1 = f"Presente/Sim (n={n1})"
    linhas = []
    for var in variaveis:
        if var not in df.columns or var == outcome_col:
            continue
        base = df[[var, outcome_col]].dropna()
        if base.empty or base[outcome_col].nunique() < 2:
            continue
        g0 = base.loc[base[outcome_col] == 0, var]
        g1 = base.loc[base[outcome_col] == 1, var]
        if len(g0) == 0 or len(g1) == 0:
            continue
        is_num = pd.api.types.is_numeric_dtype(base[var]) and base[var].nunique() > 5
        if is_num:
            g0n = pd.to_numeric(g0, errors="coerce").dropna()
            g1n = pd.to_numeric(g1, errors="coerce").dropna()
            if len(g0n) == 0 or len(g1n) == 0:
                continue
            u, p = stats.mannwhitneyu(g0n, g1n, alternative="two-sided")
            linhas.append({
                "Variável": var,
                col0: f"{g0n.median():.2f} ({g0n.quantile(0.25):.2f}–{g0n.quantile(0.75):.2f})",
                col1: f"{g1n.median():.2f} ({g1n.quantile(0.25):.2f}–{g1n.quantile(0.75):.2f})",
                "Teste": "Mann-Whitney",
                "p-valor": p,
                "Interpretação": interpretar_p(p),
            })
        else:
            tab = pd.crosstab(base[var].astype(str), base[outcome_col])
            if tab.shape[1] < 2:
                continue
            teste, p, efeito, menor_esp, chi2 = teste_categorico(tab)
            linhas.append({"Variável": var, col0: "", col1: "", "Teste": teste, "p-valor": p, "Interpretação": interpretar_p(p)})
            denom0 = max((base[outcome_col] == 0).sum(), 1)
            denom1 = max((base[outcome_col] == 1).sum(), 1)
            for nivel in tab.index:
                c0 = int(tab.loc[nivel, 0]) if 0 in tab.columns else 0
                c1 = int(tab.loc[nivel, 1]) if 1 in tab.columns else 0
                linhas.append({
                    "Variável": f"  {nivel}",
                    col0: f"{c0} ({100*c0/denom0:.1f}%)",
                    col1: f"{c1} ({100*c1/denom1:.1f}%)",
                    "Teste": "",
                    "p-valor": np.nan,
                    "Interpretação": "",
                })
    out = pd.DataFrame(linhas)
    if not out.empty:
        out["p formatado"] = out["p-valor"].apply(formatar_p)
    return out


def auditoria_qualidade(df, df_int):
    linhas = []
    n = len(df)
    def add(item, valor, denom, gravidade, interpretacao):
        linhas.append({
            "Item": item,
            "n": int(valor) if pd.notna(valor) else np.nan,
            "denominador": int(denom) if pd.notna(denom) else np.nan,
            "%": (100 * valor / denom) if denom else np.nan,
            "Gravidade": gravidade,
            "Interpretação": interpretacao,
        })
    add("Prescrições sem custo", df["CUSTO_TOTAL_R$"].isna().sum(), n, "Média", "Pode afetar análise farmacoeconômica.")
    add("DOT ausente", df["DOT_Exato"].isna().sum(), n, "Alta", "Afeta indicador de exposição terapêutica.")
    add("DOT negativo", (pd.to_numeric(df["DOT_Exato"], errors="coerce") < 0).sum(), n, "Alta", "Valor biologicamente implausível.")
    add("DDD OMS ausente/zero", ((pd.to_numeric(df["DDD_OMS"], errors="coerce").isna()) | (pd.to_numeric(df["DDD_OMS"], errors="coerce") <= 0)).sum(), n, "Média", "Impede avaliação DDD.")
    add("Outlier de dose DDD <=50% ou >=150%", df["Alerta_Outlier_Dose"].sum(), n, "Alta", "Requer revisão farmacêutica ou justificativa clínica.")
    add("AUD_PRESCRICAO não informado", df["Prescricao_Nao_Informada"].sum(), n, "Alta", "Limita inferência sobre prescrição médica.")
    add("Cultura não registrada como sim", (df["Cultura_Sim"] == 0).sum(), n, "Média", "Limita análise microbiológica/descalonamento.")
    add("Tempo de tratamento não determinado", (df["Tempo_Determinado"] == 0).sum(), n, "Média", "Afeta análise de conformidade terapêutica.")
    add("Internações com múltiplos desfechos", df_int["multiplos_desfechos"].sum(), len(df_int), "Alta", "Requer regra de hierarquia ou validação manual.")
    return pd.DataFrame(linhas)


def calculate_diagnostic_metrics(sensitivity, specificity, prevalence):
    eps = 1e-12
    ppv = (sensitivity * prevalence) / ((sensitivity * prevalence) + (1 - specificity) * (1 - prevalence) + eps)
    npv = (specificity * (1 - prevalence)) / ((specificity * (1 - prevalence)) + (1 - sensitivity) * prevalence + eps)
    return ppv, npv

# ============================================================
# 5. LEITURA DOS DADOS
# ============================================================

uploaded = st.sidebar.file_uploader("Enviar planilha Excel ou CSV", type=["xlsx", "xls", "csv"])

try:
    df_original = carregar_arquivo(uploaded)
    df = preparar_variaveis(df_original)
    df_internacao = agregar_internacoes(df)
    df_sobreposicao = detectar_sobreposicao_potencial(df)
except Exception as e:
    st.error(f"Erro ao carregar dados: {e}")
    st.stop()

# ============================================================
# 6. FILTROS
# ============================================================

st.sidebar.header("Filtros")
setores = sorted(df["Setor_Padronizado"].dropna().astype(str).unique().tolist())
setor_selecionado = st.sidebar.multiselect("Setor de internação", options=setores, default=setores)

motivos = sorted(df["Motivo"].dropna().astype(str).unique().tolist())
motivo_selecionado = st.sidebar.multiselect("Motivo da prescrição", options=motivos, default=motivos)

aware_opts = sorted(df["AWaRe"].dropna().astype(str).unique().tolist())
aware_selecionado = st.sidebar.multiselect("AWaRe", options=aware_opts, default=aware_opts)

df_filtrado = df[
    df["Setor_Padronizado"].astype(str).isin(setor_selecionado)
    & df["Motivo"].astype(str).isin(motivo_selecionado)
    & df["AWaRe"].astype(str).isin(aware_selecionado)
].copy()

ids_filtrados = df_filtrado["ID_Internacao"].dropna().unique().tolist()
df_internacao_filtrado = df_internacao[df_internacao["ID_Internacao"].isin(ids_filtrados)].copy()
df_sobreposicao_filtrado = df_sobreposicao[df_sobreposicao["ID_Internacao"].isin(ids_filtrados)].copy() if not df_sobreposicao.empty else df_sobreposicao

st.caption(
    "Convenção estatística: p<0,05 sugere evidência de associação/diferença; "
    "p≥0,05 indica ausência de evidência suficiente, não prova ausência de efeito. "
    "As análises de mortalidade e cenários são exploratórias quando a amostra é pequena."
)

# ============================================================
# 7. ABAS DO DASHBOARD
# ============================================================

tab1, tab2, tab3, tab4, tab5, tab6, tab7, tab8, tab9, tab10 = st.tabs([
    "👁️ 1. Visão Geral",
    "📈 2. Padrões Clínicos",
    "🏥 3. Desfechos e Risco",
    "💰 4. Farmacoeconomia",
    "🎯 5. Auditoria",
    "🧬 6. Stewardship Avançado",
    "🗄️ 7. Dados",
    "📋 8. Qualidade & Artigo",
    "🔮 9. Modelagem & Bootstrap",
    "🧪 10. Sensibilidade",
])

# ------------------------------------------------------------
# ABA 1
# ------------------------------------------------------------
with tab1:
    c1, c2, c3, c4, c5, c6 = st.columns(6)
    c1.metric("Internações", f"{df_internacao_filtrado['ID_Internacao'].nunique()}")
    c2.metric("Prescrições", f"{len(df_filtrado)}")
    c3.metric("Antimicrobianos", f"{df_filtrado['ATB'].nunique()}")
    c4.metric("Custo total", f"R$ {df_filtrado['CUSTO_TOTAL_R$'].sum():,.2f}")
    c5.metric("DOT mediano", f"{df_filtrado['DOT_Exato'].median():.2f}")
    c6.metric("ASI médio", f"{df_filtrado['ASI'].mean():.2f}" if df_filtrado['ASI'].notna().any() else "NA")

    st.markdown("---")
    col1, col2 = st.columns(2)
    with col1:
        st.plotly_chart(px.pie(df_filtrado, names="Setor_Padronizado", title="Prescrições por setor", hole=0.45), use_container_width=True)
    with col2:
        top_atb = df_filtrado["ATB"].value_counts().reset_index()
        top_atb.columns = ["ATB", "n"]
        st.plotly_chart(px.bar(top_atb.head(15), x="n", y="ATB", orientation="h", title="Antimicrobianos mais prescritos"), use_container_width=True)

    col3, col4 = st.columns(2)
    with col3:
        st.plotly_chart(px.pie(df_filtrado, names="Motivo", title="Motivo da prescrição", hole=0.45), use_container_width=True)
    with col4:
        st.plotly_chart(px.pie(df_filtrado, names="AWaRe", title="Classificação AWaRe aproximada", hole=0.45), use_container_width=True)

# ------------------------------------------------------------
# ABA 2
# ------------------------------------------------------------
with tab2:
    st.subheader("Padrões clínicos: DDD, DOT, classes, síndromes e ASI")
    col1, col2 = st.columns(2)
    with col1:
        st.plotly_chart(px.pie(df_filtrado, names="Status_DDD", title="Classificação pela DDD", hole=0.45), use_container_width=True)
    with col2:
        st.plotly_chart(px.box(df_filtrado, x="Setor_Padronizado", y="DOT_Exato", points="all", title="DOT por setor"), use_container_width=True)

    testar_kruskal_e_posthoc(df_filtrado, "Setor_Padronizado", "DOT_Exato", "Kruskal-Wallis para DOT por setor")
    testar_kruskal_e_posthoc(df_filtrado.dropna(subset=["ASI"]), "Setor_Padronizado", "ASI", "Kruskal-Wallis para ASI por setor")

    st.markdown("---")
    col3, col4 = st.columns(2)
    with col3:
        tabela_classes = pd.crosstab(df_filtrado["Setor_Padronizado"], df_filtrado["CLASSE"])
        st.plotly_chart(px.bar(tabela_classes.reset_index(), x="Setor_Padronizado", y=tabela_classes.columns, barmode="stack", title="Classes farmacológicas por setor"), use_container_width=True)
        teste, p, efeito, menor_esp, chi2 = teste_categorico(tabela_classes)
        detalhes = f"efeito={efeito:.2f}".replace(".", ",") if pd.notna(efeito) else ""
        if pd.notna(menor_esp):
            detalhes += f"; menor esperado={menor_esp:.2f}".replace(".", ",")
        mostrar_resultado_teste(f"{teste} para classe farmacológica por setor", estatistica=chi2, p=p, detalhes=detalhes)
    with col4:
        top_sind = df_filtrado["Sindrome"].value_counts().reset_index()
        top_sind.columns = ["Síndrome", "n"]
        st.plotly_chart(px.bar(top_sind, x="n", y="Síndrome", orientation="h", title="Síndromes/indicações"), use_container_width=True)
        tabela_ddd_setor = pd.crosstab(df_filtrado["Setor_Padronizado"], df_filtrado["Status_DDD"])
        st.write("Status DDD por setor")
        st.dataframe(tabela_ddd_setor, use_container_width=True)

# ------------------------------------------------------------
# ABA 3
# ------------------------------------------------------------
with tab3:
    st.subheader("Desfechos e risco: análise por internação")
    st.caption("Esta aba usa base agregada por internação para reduzir pseudorrepetição.")
    col1, col2 = st.columns(2)
    with col1:
        tabela_desf = pd.crosstab(df_internacao_filtrado["setor_principal"], df_internacao_filtrado["desfecho_final"])
        st.plotly_chart(px.bar(tabela_desf.reset_index(), x="setor_principal", y=tabela_desf.columns, barmode="group", title="Desfecho final por setor principal"), use_container_width=True)
        teste, p, efeito, menor_esp, chi2 = teste_categorico(tabela_desf)
        detalhes = f"efeito={efeito:.2f}".replace(".", ",") if pd.notna(efeito) else ""
        if pd.notna(menor_esp):
            detalhes += f"; menor esperado={menor_esp:.2f}".replace(".", ",")
        mostrar_resultado_teste(f"{teste} para desfecho final por setor", estatistica=chi2, p=p, detalhes=detalhes)
    with col2:
        cols_resumo = ["ID_Internacao", "setor_principal", "idade", "n_prescricoes", "n_atb", "dot_total", "custo_total", "asi_total", "desfecho_final"]
        st.dataframe(df_internacao_filtrado[cols_resumo], use_container_width=True)

    st.markdown("---")
    st.subheader("Odds Ratio exploratório 2x2")
    col3, col4, col5 = st.columns(3)
    with col3:
        exposicao = st.selectbox("Exposição", ["amplo_espectro", "erro_admin_algum", "tempo_determinado_algum", "cultura_alguma", "aware_watch_reserve", "prof_maior_24h_alguma", "outlier_dose_algum"], index=0)
    with col4:
        desfecho_bin = st.selectbox("Desfecho binário", ["obito", "erro_admin_algum", "amplo_espectro", "prof_maior_24h_alguma"], index=0)
    with col5:
        st.caption("OR bruto com Fisher exato. Não interpretar como causalidade.")
    if exposicao == desfecho_bin:
        st.warning("Escolha exposição e desfecho diferentes.")
    else:
        tab_or = tabela_2x2_binaria(df_internacao_filtrado, exposicao, desfecho_bin)
        st.dataframe(tab_or, use_container_width=True)
        odds, li, ls, p = fisher_or_com_ic(tab_or)
        detalhes = f"OR={odds:.2f}; IC95% {li:.2f}–{ls:.2f}".replace(".", ",")
        mostrar_resultado_teste("Fisher exato para OR bruto", p=p, detalhes=detalhes)

# ------------------------------------------------------------
# ABA 4
# ------------------------------------------------------------
with tab4:
    st.subheader("Farmacoeconomia: custo observado e custo associado a inconformidades")
    col1, col2 = st.columns(2)
    with col1:
        st.plotly_chart(px.box(df_filtrado, x="Setor_Padronizado", y="CUSTO_TOTAL_R$", points="all", title="Custo por prescrição"), use_container_width=True)
    with col2:
        st.plotly_chart(px.scatter(df_filtrado, x="DOT_Exato", y="CUSTO_TOTAL_R$", color="Setor_Padronizado", trendline="ols", title="Custo por DOT - prescrição"), use_container_width=True)
    testar_kruskal_e_posthoc(df_filtrado, "Setor_Padronizado", "CUSTO_TOTAL_R$", "Kruskal-Wallis para custo por setor")

    st.markdown("---")
    k1, k2, k3 = st.columns(3)
    custo_total = df_filtrado["CUSTO_TOTAL_R$"].sum()
    custo_inconf = df_filtrado["Custo_Inconformidade"].sum()
    k1.metric("Custo total", f"R$ {custo_total:,.2f}")
    k2.metric("Custo associado a inconformidades", f"R$ {custo_inconf:,.2f}")
    k3.metric("Percentual sensível à intervenção", f"{100*custo_inconf/custo_total:.1f}%" if custo_total > 0 else "NA")
    st.caption("Não interpretar como economia real; trata-se de custo associado a inconformidades potencialmente sensíveis à intervenção de stewardship.")

    col3, col4 = st.columns(2)
    with col3:
        st.plotly_chart(px.box(df_internacao_filtrado, x="setor_principal", y="custo_total", points="all", title="Custo total por internação"), use_container_width=True)
    with col4:
        st.plotly_chart(px.scatter(df_internacao_filtrado, x="dot_total", y="custo_total", color="setor_principal", size="n_atb", title="Custo total x DOT total por internação"), use_container_width=True)
    corr_custo = correlacoes_spearman(df_internacao_filtrado, ["dot_total", "custo_total", "n_atb", "n_prescricoes", "asi_total"])
    st.dataframe(corr_custo.round(4), use_container_width=True)

# ------------------------------------------------------------
# ABA 5
# ------------------------------------------------------------
with tab5:
    st.subheader("Auditoria: associações categóricas e falhas operacionais")
    col1, col2 = st.columns(2)
    with col1:
        tb_admin_setor = pd.crosstab(df_filtrado["Setor_Padronizado"], df_filtrado["Erro_Administracao"])
        st.write("Erro/incompletude de administração por setor")
        st.dataframe(tb_admin_setor, use_container_width=True)
        teste, p, efeito, menor_esp, chi2 = teste_categorico(tb_admin_setor)
        mostrar_resultado_teste(f"{teste}: erro de administração por setor", estatistica=chi2, p=p, detalhes=f"efeito={efeito:.2f}".replace(".", ",") if pd.notna(efeito) else "")
    with col2:
        tb_motivo = pd.crosstab(df_filtrado["Motivo"], df_filtrado["Erro_Administracao"])
        st.write("Erro/incompletude de administração por motivo")
        st.dataframe(tb_motivo, use_container_width=True)
        teste, p, efeito, menor_esp, chi2 = teste_categorico(tb_motivo)
        mostrar_resultado_teste(f"{teste}: erro de administração por motivo", estatistica=chi2, p=p, detalhes=f"efeito={efeito:.2f}".replace(".", ",") if pd.notna(efeito) else "")

    st.markdown("---")
    col3, col4 = st.columns(2)
    with col3:
        tb_ddd = pd.crosstab(df_filtrado["Status_DDD"], df_filtrado["Erro_Administracao"])
        st.write("Status DDD x erro de administração")
        st.dataframe(tb_ddd, use_container_width=True)
        teste, p, efeito, menor_esp, chi2 = teste_categorico(tb_ddd)
        mostrar_resultado_teste(f"{teste}: DDD x erro de administração", estatistica=chi2, p=p)
    with col4:
        tb_aware = pd.crosstab(df_filtrado["AWaRe"], df_filtrado["Erro_Administracao"])
        st.write("AWaRe x erro de administração")
        st.dataframe(tb_aware, use_container_width=True)
        teste, p, efeito, menor_esp, chi2 = teste_categorico(tb_aware)
        mostrar_resultado_teste(f"{teste}: AWaRe x erro de administração", estatistica=chi2, p=p)

# ------------------------------------------------------------
# ABA 6 - STEWARDSHIP AVANÇADO
# ------------------------------------------------------------
with tab6:
    st.subheader("Stewardship avançado sem alterar a planilha")
    st.info("As análises abaixo são derivadas dos campos já existentes. Revisar clinicamente antes de usar como resultado definitivo de artigo.")

    st.markdown("### 1. Ranking de antimicrobianos por impacto")
    ranking_atb = df_filtrado.groupby("ATB", dropna=False).agg(
        n=("ATB", "size"),
        dot_total=("DOT_Exato", "sum"),
        ddd_total=("DDD_TOTAL", "sum"),
        custo_total=("CUSTO_TOTAL_R$", "sum"),
        erro_admin_pct=("Erro_Administracao", "mean"),
        ddd_inadequada_pct=("DDD_Inadequada", "mean"),
        asi_medio=("ASI", "mean"),
    ).reset_index()
    ranking_atb["erro_admin_pct"] *= 100
    ranking_atb["ddd_inadequada_pct"] *= 100
    ranking_atb = ranking_atb.sort_values(["custo_total", "dot_total", "n"], ascending=False)
    st.dataframe(ranking_atb.round(3), use_container_width=True)
    st.download_button("Baixar ranking de antimicrobianos", ranking_atb.to_csv(index=False, sep=";", decimal=",").encode("utf-8-sig"), "ranking_antimicrobianos.csv", "text/csv")

    st.markdown("---")
    st.markdown("### 2. AWaRe e Índice de Espectro Antimicrobiano (ASI)")
    col1, col2 = st.columns(2)
    with col1:
        aware_tab = pd.crosstab(df_filtrado["Setor_Padronizado"], df_filtrado["AWaRe"], normalize="index") * 100
        st.plotly_chart(px.bar(aware_tab.reset_index(), x="Setor_Padronizado", y=aware_tab.columns, barmode="stack", title="AWaRe por setor (%)"), use_container_width=True)
    with col2:
        st.plotly_chart(px.box(df_filtrado.dropna(subset=["ASI"]), x="Setor_Padronizado", y="ASI", points="all", title="ASI por setor"), use_container_width=True)

    st.markdown("---")
    st.markdown("### 3. Profilaxia cirúrgica")
    prof = df_filtrado[df_filtrado["Profilaxia_Cirurgica"] == 1].copy()
    if prof.empty:
        st.warning("Não foram identificadas prescrições classificáveis como profilaxia cirúrgica pelos campos existentes.")
    else:
        k1, k2, k3 = st.columns(3)
        k1.metric("Profilaxias identificadas", len(prof))
        k2.metric(">24h", f"{prof['Profilaxia_Maior_24h'].sum()} ({100*prof['Profilaxia_Maior_24h'].mean():.1f}%)")
        k3.metric(">48h", f"{prof['Profilaxia_Maior_48h'].sum()} ({100*prof['Profilaxia_Maior_48h'].mean():.1f}%)")
        st.plotly_chart(px.box(prof, x="Setor_Padronizado", y="Tempo_Analise_Profilaxia", points="all", title="Duração estimada da profilaxia cirúrgica"), use_container_width=True)
        tb_prof = pd.crosstab(prof["Setor_Padronizado"], prof["Profilaxia_Maior_24h"])
        st.write("Profilaxia >24h por setor")
        st.dataframe(tb_prof, use_container_width=True)

    st.markdown("---")
    st.markdown("### 4. Double coverage / sobreposição terapêutica potencial")
    st.caption("Quando há datas de início/fim, o app verifica sobreposição temporal. Sem datas confiáveis, sinaliza apenas sobreposição potencial por mesma classe/espectro na internação.")
    st.dataframe(df_sobreposicao_filtrado, use_container_width=True)
    st.download_button("Baixar possíveis sobreposições", df_sobreposicao_filtrado.to_csv(index=False, sep=";", decimal=",").encode("utf-8-sig"), "sobreposicao_potencial.csv", "text/csv")

# ------------------------------------------------------------
# ABA 7 - DADOS
# ------------------------------------------------------------
with tab7:
    st.subheader("Base de prescrições filtrada")
    st.dataframe(df_filtrado, use_container_width=True, height=420)
    st.download_button("Baixar prescrições tratadas em CSV", df_filtrado.to_csv(index=False, sep=";", decimal=",").encode("utf-8-sig"), "prescricoes_tratadas.csv", "text/csv")
    st.subheader("Base agregada por internação")
    st.dataframe(df_internacao_filtrado, use_container_width=True, height=320)
    st.download_button("Baixar base por internação em CSV", df_internacao_filtrado.to_csv(index=False, sep=";", decimal=",").encode("utf-8-sig"), "internacoes_agregadas.csv", "text/csv")

# ------------------------------------------------------------
# ABA 8 - QUALIDADE E ARTIGO
# ------------------------------------------------------------
with tab8:
    st.subheader("Qualidade do banco, Tabela 1 e saídas para artigo")
    n_presc = len(df_filtrado)
    n_int = df_internacao_filtrado["ID_Internacao"].nunique()
    n_obitos = int(df_internacao_filtrado["obito"].sum()) if len(df_internacao_filtrado) else 0
    n_admin_erro = int(df_filtrado["Erro_Administracao"].sum()) if n_presc else 0
    n_ddd_ok = int(df_filtrado["DDD_Adequada"].sum()) if n_presc else 0
    n_aware_wr = int(df_filtrado["AWaRe"].isin(["Watch", "Reserve"]).sum()) if n_presc else 0

    k1, k2, k3, k4, k5 = st.columns(5)
    k1.metric("Prescrições", n_presc)
    k2.metric("Internações", n_int)
    k3.metric("Óbitos", n_obitos)
    k4.metric("Erro admin.", formatar_pct_ic(n_admin_erro, n_presc))
    k5.metric("Watch/Reserve", formatar_pct_ic(n_aware_wr, n_presc))

    indicadores = pd.DataFrame([
        ["Adequação à DDD", n_ddd_ok, n_presc, formatar_pct_ic(n_ddd_ok, n_presc)],
        ["Administração incompleta/suspensa", n_admin_erro, n_presc, formatar_pct_ic(n_admin_erro, n_presc)],
        ["Tempo determinado pelo prescritor", int(df_filtrado["Tempo_Determinado"].sum()), n_presc, formatar_pct_ic(int(df_filtrado["Tempo_Determinado"].sum()), n_presc)],
        ["Cultura registrada como sim", int(df_filtrado["Cultura_Sim"].sum()), n_presc, formatar_pct_ic(int(df_filtrado["Cultura_Sim"].sum()), n_presc)],
        ["AUD_PRESCRICAO não informado", int(df_filtrado["Prescricao_Nao_Informada"].sum()), n_presc, formatar_pct_ic(int(df_filtrado["Prescricao_Nao_Informada"].sum()), n_presc)],
        ["Óbito por internação", n_obitos, n_int, formatar_pct_ic(n_obitos, n_int)],
        ["Outlier de dose", int(df_filtrado["Alerta_Outlier_Dose"].sum()), n_presc, formatar_pct_ic(int(df_filtrado["Alerta_Outlier_Dose"].sum()), n_presc)],
    ], columns=["Indicador", "n", "denominador", "% e IC95%"])
    st.dataframe(indicadores, use_container_width=True)

    st.markdown("---")
    st.markdown("### Auditoria de qualidade dos dados")
    qualidade = auditoria_qualidade(df_filtrado, df_internacao_filtrado)
    st.dataframe(qualidade.round(2), use_container_width=True)

    st.markdown("---")
    st.markdown("### Tabela 1 automática")
    outcome = st.selectbox("Estratificar Tabela 1 por", ["erro_admin_algum", "obito", "amplo_espectro", "prof_maior_24h_alguma"], index=0)
    vars_t1 = ["idade", "sexo", "setor_principal", "n_prescricoes", "n_atb", "dot_total", "ddd_total", "custo_total", "asi_total", "amplo_espectro", "aware_watch_reserve", "tempo_determinado_algum", "cultura_alguma", "outlier_dose_algum"]
    tabela1 = gerar_tabela1(df_internacao_filtrado, outcome, vars_t1)
    if tabela1.empty:
        st.warning("Não foi possível gerar Tabela 1 para o desfecho selecionado.")
    else:
        st.dataframe(tabela1, use_container_width=True)
        st.download_button("Baixar Tabela 1", tabela1.to_csv(index=False, sep=";", decimal=",").encode("utf-8-sig"), "tabela1_artigo.csv", "text/csv")

    st.markdown("---")
    st.write("Texto-base para método/resultados")
    st.code(
        f"""Estudo observacional retrospectivo baseado em {n_presc} prescrições antimicrobianas referentes a {n_int} internações.
As análises foram conduzidas em duas unidades: prescrição e internação. Para variáveis contínuas assimétricas, utilizaram-se mediana e intervalo interquartil, com comparação por Kruskal-Wallis ou Mann-Whitney. Proporções foram apresentadas com IC95% pelo método de Wilson. Variáveis categóricas foram avaliadas por Fisher exato em tabelas 2x2 ou qui-quadrado quando aplicável.
A adequação à DDD foi observada em {formatar_pct_ic(n_ddd_ok, n_presc)} das prescrições. Administração incompleta ou suspensão precoce ocorreu em {formatar_pct_ic(n_admin_erro, n_presc)}. Antimicrobianos Watch/Reserve representaram {formatar_pct_ic(n_aware_wr, n_presc)}. Os resultados devem ser interpretados como associações, sem inferência causal.""",
        language="markdown",
    )

# ------------------------------------------------------------
# ABA 9 - MODELAGEM & BOOTSTRAP
# ------------------------------------------------------------
with tab9:
    st.subheader("Modelagem, correlação, bootstrap e limites inferenciais")
    st.markdown("### 1. Correlação de Spearman")
    vars_corr_presc = ["Idade", "DOT_Exato", "CUSTO_TOTAL_R$", "Ratio_DDD", "ASI"]
    interp_corr = correlacoes_spearman(df_filtrado, vars_corr_presc)
    st.dataframe(interp_corr.round(4), use_container_width=True)

    st.markdown("---")
    st.markdown("### 2. Bootstrap para medianas")
    variavel_boot = st.selectbox("Variável para IC95% bootstrap da mediana", ["DOT_Exato", "CUSTO_TOTAL_R$", "Ratio_DDD", "ASI"])
    med, li, ls = bootstrap_ci_median(df_filtrado[variavel_boot])
    st.metric("Mediana observada", f"{med:.3f}" if pd.notna(med) else "NA")
    st.info(f"IC95% bootstrap da mediana: {li:.3f} a {ls:.3f}" if pd.notna(li) else "Bootstrap não estimável.")

    st.markdown("### 3. Bootstrap para Spearman")
    colx, coly = st.columns(2)
    with colx:
        xboot = st.selectbox("X", ["DOT_Exato", "CUSTO_TOTAL_R$", "Ratio_DDD", "ASI", "Idade"], index=0)
    with coly:
        yboot = st.selectbox("Y", ["CUSTO_TOTAL_R$", "DOT_Exato", "Ratio_DDD", "ASI", "Idade"], index=0)
    rho, rli, rls = bootstrap_ci_spearman(df_filtrado, xboot, yboot)
    st.info(f"rho={rho:.3f}; IC95% bootstrap {rli:.3f} a {rls:.3f}" if pd.notna(rho) else "Correlação não estimável.")

    st.markdown("---")
    st.markdown("### 4. Mortalidade — análise exploratória")
    eventos = int(df_internacao_filtrado["obito"].sum()) if len(df_internacao_filtrado) else 0
    n_modelo = len(df_internacao_filtrado)
    st.write(f"Eventos de óbito por internação: **{eventos}/{n_modelo}**.")
    if eventos < 10:
        st.warning("Regressão logística múltipla não é recomendada com menos de 10 eventos. A análise abaixo é exploratória.")
        comparacoes = []
        for var in ["idade", "dot_total", "custo_total", "n_prescricoes", "n_atb", "asi_total", "custo_inconformidade"]:
            g0 = pd.to_numeric(df_internacao_filtrado.loc[df_internacao_filtrado["obito"] == 0, var], errors="coerce").dropna()
            g1 = pd.to_numeric(df_internacao_filtrado.loc[df_internacao_filtrado["obito"] == 1, var], errors="coerce").dropna()
            if len(g0) > 0 and len(g1) > 0:
                u, p = stats.mannwhitneyu(g1, g0, alternative="two-sided")
                comparacoes.append({"Variável": var, "Mediana óbito": g1.median(), "Mediana não óbito": g0.median(), "U": u, "p-valor": p, "p formatado": formatar_p(p), "Interpretação": interpretar_p(p)})
        st.dataframe(pd.DataFrame(comparacoes).round(4), use_container_width=True)
    else:
        base = df_internacao_filtrado[["obito", "idade", "amplo_espectro", "erro_admin_algum", "n_atb", "asi_total"]].dropna()
        X = sm.add_constant(base[["idade", "amplo_espectro", "erro_admin_algum", "n_atb", "asi_total"]])
        y = base["obito"]
        try:
            modelo = sm.Logit(y, X).fit(disp=False)
            res = pd.DataFrame({"OR ajustado": np.exp(modelo.params), "IC95% inferior": np.exp(modelo.conf_int()[0]), "IC95% superior": np.exp(modelo.conf_int()[1]), "p-valor": modelo.pvalues})
            res["p formatado"] = res["p-valor"].apply(formatar_p)
            res["Interpretação"] = res["p-valor"].apply(interpretar_p)
            st.dataframe(res.drop("const", errors="ignore"), use_container_width=True)
        except Exception as e:
            st.warning(f"O modelo não convergiu: {e}")
    st.caption("Firth/Pydantic não foram ativados nesta versão para evitar novas dependências e manter publicação simples no Streamlit. A análise de eventos raros permanece protegida por trava metodológica.")

# ------------------------------------------------------------
# ABA 10 - SENSIBILIDADE
# ------------------------------------------------------------
with tab10:
    st.subheader("Análise de sensibilidade epidemiológica — cenários simulados")
    st.warning("Esta aba é uma simulação didática/estratégica. Não usa dados reais de resistência microbiológica da planilha atual e não deve ser descrita como resultado clínico.")
    col1, col2, col3 = st.columns(3)
    with col1:
        sens = st.slider("Sensibilidade do alerta/modelo", 0.01, 0.99, 0.85, 0.01)
    with col2:
        spec = st.slider("Especificidade do alerta/modelo", 0.01, 0.99, 0.90, 0.01)
    with col3:
        prev_atual = st.slider("Prevalência simulada de MDRO", 0.01, 0.80, 0.15, 0.01)
    ppv, npv = calculate_diagnostic_metrics(sens, spec, prev_atual)
    m1, m2 = st.columns(2)
    m1.metric("VPP simulado", f"{100*ppv:.1f}%")
    m2.metric("VPN simulado", f"{100*npv:.1f}%")

    prevalencias = np.linspace(0.01, 0.80, 80)
    df_sim = pd.DataFrame([
        {"Prevalência": p, "VPP": calculate_diagnostic_metrics(sens, spec, p)[0], "VPN": calculate_diagnostic_metrics(sens, spec, p)[1]}
        for p in prevalencias
    ])
    st.plotly_chart(px.line(df_sim, x="Prevalência", y=["VPP", "VPN"], title="Variação de VPP/VPN conforme prevalência simulada"), use_container_width=True)
    st.code(
        """Uso recomendado na dissertação:
A análise de sensibilidade epidemiológica foi incluída como funcionalidade exploratória do dashboard para simular o impacto da prevalência de patógenos multirresistentes sobre valores preditivos de regras de alerta. Como a base atual não contém dados microbiológicos estruturados suficientes, essa aba não foi utilizada como inferência clínica principal.""",
        language="markdown",
    )
