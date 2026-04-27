# -*- coding: utf-8 -*-
"""
Dashboard Executivo & Estatístico - Stewardship Antimicrobiano
Mestrado Fiocruz - Versão completa V17

Como usar:
1) Coloque este arquivo .py na mesma pasta da planilha:
   Sistema de Monitoramento - Mestrado Jobson (V14.0).xlsx
2) Execute:
   python -m streamlit run app_mestrado_jobson_v17_completo.py

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
    page_title="Stewardship Antimicrobiano - Mestrado",
    page_icon="💊",
    layout="wide"
)

st.title("📊 Dashboard Executivo & Estatístico - Stewardship Antimicrobiano")
st.caption(
    "Versão completa para dissertação/artigo: leitura da aba PRESCRIÇÕES, "
    "análise por prescrição e por internação, IC95%, testes não paramétricos, "
    "associação categórica, correlação, pós-testes e interpretações automáticas dos p-values."
)
st.markdown("---")

ARQUIVO_PADRAO = "Sistema de Monitoramento - Mestrado Jobson (V14.0).xlsx"
ABA_PADRAO = "PRESCRIÇÕES"
ALFA = 0.05

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

# ============================================================
# 2. FUNÇÕES DE UTILIDADE TEXTUAL E ESTATÍSTICA
# ============================================================


def remover_acentos(texto):
    """Remove acentos e normaliza textos para busca/padronização."""
    if pd.isna(texto):
        return ""
    texto = str(texto).strip()
    return "".join(
        ch for ch in unicodedata.normalize("NFKD", texto)
        if not unicodedata.combining(ch)
    )


def padronizar_setor(valor):
    """Agrupa variações de escrita dos setores."""
    x = remover_acentos(valor).upper()
    x = re.sub(r"[\-_/]+", " ", x)
    x = re.sub(r"\s+", " ", x).strip()

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
    """Converte datas serial do Excel para datetime."""
    return pd.to_datetime(serie, unit="D", origin="1899-12-30", errors="coerce")


def classificar_ddd(ratio):
    """Classifica dose diária prescrita em relação à DDD OMS."""
    if pd.isna(ratio) or ratio <= 0:
        return "Não Avaliado"
    if 0.9 <= ratio <= 1.1:
        return "Adequada à DDD (±10%)"
    if ratio < 0.9:
        return "Subdose"
    return "Superdose"


def wilson_ic(x, n, alpha=0.05):
    """Intervalo de confiança de Wilson para proporções."""
    if n == 0 or pd.isna(n):
        return np.nan, np.nan
    z = stats.norm.ppf(1 - alpha / 2)
    p = x / n
    den = 1 + z**2 / n
    centro = (p + z**2 / (2 * n)) / den
    margem = z * np.sqrt((p * (1 - p) + z**2 / (4 * n)) / n) / den
    return centro - margem, centro + margem


def formatar_pct_ic(x, n):
    """Formata proporção com IC95% de Wilson."""
    if n == 0 or pd.isna(n):
        return "NA"
    li, ls = wilson_ic(x, n)
    return f"{100*x/n:.1f}% (IC95% {100*li:.1f}–{100*ls:.1f})"


def formatar_p(p):
    """Formata p-valor em padrão brasileiro."""
    if p is None or pd.isna(p):
        return "NA"
    if p < 0.001:
        return "<0,001"
    return f"{p:.4f}".replace(".", ",")


def interpretar_p(p, alfa=ALFA):
    """Interpreta p-valor sem sugerir causalidade."""
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
    """Classifica força e direção da correlação de Spearman."""
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
    """Mostra teste estatístico com p-valor e interpretação padronizada."""
    partes = [f"**{nome}**"]
    if estatistica is not None and pd.notna(estatistica):
        partes.append(f"estatística={estatistica:.3f}".replace(".", ","))
    if p is not None and pd.notna(p):
        partes.append(f"p-valor={formatar_p(p)}")
    if detalhes:
        partes.append(detalhes)

    texto = "; ".join(partes) + "."
    st.info(texto)

    interp = interpretar_p(p, alfa)
    if p is not None and pd.notna(p) and p < alfa:
        st.success(interp)
    else:
        st.warning(interp)


def descrever_num(serie):
    """Resumo robusto para variáveis numéricas."""
    s = pd.to_numeric(serie, errors="coerce").dropna()
    if len(s) == 0:
        return pd.Series({
            "n": 0, "média": np.nan, "DP": np.nan, "mediana": np.nan,
            "Q1": np.nan, "Q3": np.nan, "mín": np.nan, "máx": np.nan
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
    """Define desfecho final por hierarquia clínica quando há múltiplos registros."""
    vals = grupo.dropna().astype(str).unique().tolist()
    if not vals:
        return np.nan
    vals_ordenados = sorted(vals, key=lambda v: HIERARQUIA_DESFECHO.get(v, 0), reverse=True)
    return vals_ordenados[0]


def moda_segura(serie):
    m = serie.dropna().mode()
    return m.iloc[0] if len(m) else np.nan


def nome_coluna_prescricao(df):
    """Retorna o nome da coluna de prescrição, com ou sem acento."""
    if "ID_PRESCRIÇÃO" in df.columns:
        return "ID_PRESCRIÇÃO"
    if "ID_PRESCRICAO" in df.columns:
        return "ID_PRESCRICAO"
    return None

# ============================================================
# 3. CARREGAMENTO E PREPARAÇÃO DO BANCO
# ============================================================


def carregar_arquivo(uploaded_file=None):
    """Carrega Excel ou CSV. Por padrão, usa a aba PRESCRIÇÕES do Excel."""
    if uploaded_file is not None:
        nome = uploaded_file.name.lower()
        if nome.endswith((".xlsx", ".xls")):
            xls = pd.ExcelFile(uploaded_file)
            aba = ABA_PADRAO if ABA_PADRAO in xls.sheet_names else xls.sheet_names[0]
            return pd.read_excel(uploaded_file, sheet_name=aba)
        return pd.read_csv(uploaded_file, sep=";", decimal=",")

    caminho_excel = Path(ARQUIVO_PADRAO)
    if caminho_excel.exists():
        xls = pd.ExcelFile(caminho_excel)
        aba = ABA_PADRAO if ABA_PADRAO in xls.sheet_names else xls.sheet_names[0]
        return pd.read_excel(caminho_excel, sheet_name=aba)

    caminho_csv = Path("PRESCRIÇÕES.csv")
    if caminho_csv.exists():
        return pd.read_csv(caminho_csv, sep=";", decimal=",")

    raise FileNotFoundError(
        f"Coloque o arquivo '{ARQUIVO_PADRAO}' na mesma pasta do app, "
        "ou envie a planilha pela barra lateral."
    )


def preparar_variaveis(df):
    """Limpa, padroniza e cria variáveis analíticas."""
    df = df.copy()
    df.columns = [str(c).strip() for c in df.columns]

    if "ID_Internacao" not in df.columns:
        raise ValueError("A coluna 'ID_Internacao' não foi encontrada.")

    df = df.dropna(subset=["ID_Internacao"]).copy()

    for col in COLUNAS_NUMERICAS:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")

    for col in ["Carimbo", "Inicio", "Fim"]:
        if col in df.columns:
            numerica = pd.to_numeric(df[col], errors="coerce")
            if numerica.notna().sum() > 0:
                df[col + "_Data"] = excel_serial_para_data(numerica)

    if "Setor" in df.columns:
        df["Setor_Padronizado"] = df["Setor"].apply(padronizar_setor)
    else:
        df["Setor_Padronizado"] = np.nan

    if "Idade" in df.columns:
        bins = [0, 17, 39, 59, 130]
        labels = [
            "1. Pediátrico (<18)",
            "2. Adulto jovem (18–39)",
            "3. Adulto (40–59)",
            "4. Idoso (60+)",
        ]
        df["Faixa_Etaria"] = pd.cut(df["Idade"], bins=bins, labels=labels, right=True)

    if {"Dose_diaria", "DDD_OMS"}.issubset(df.columns):
        df["Ratio_DDD"] = np.where(
            df["DDD_OMS"] > 0,
            (df["Dose_diaria"] / 1000) / df["DDD_OMS"],
            np.nan,
        )
        df["Status_DDD"] = df["Ratio_DDD"].apply(classificar_ddd)
        df["DDD_Adequada"] = df["Status_DDD"].eq("Adequada à DDD (±10%)").astype(int)
        df["DDD_Inadequada"] = df["Status_DDD"].isin(["Subdose", "Superdose"]).astype(int)
    else:
        df["Ratio_DDD"] = np.nan
        df["Status_DDD"] = "Não Avaliado"
        df["DDD_Adequada"] = 0
        df["DDD_Inadequada"] = 0

    aud_presc = df.get("AUD_PRESCRICAO", pd.Series([""] * len(df), index=df.index)).astype(str)
    aud_admin = df.get("AUD_ADMINISTRACAO", pd.Series([""] * len(df), index=df.index)).astype(str)
    desfecho = df.get("Desfecho", pd.Series([""] * len(df), index=df.index)).astype(str)
    cultura = df.get("Cultura", pd.Series([""] * len(df), index=df.index)).astype(str)
    espectro = df.get("ESPECTRO", pd.Series([""] * len(df), index=df.index)).astype(str)
    tempo_det = df.get("Determinou_Tempo?", pd.Series([""] * len(df), index=df.index)).astype(str)

    # Prescrição: não tratar "Não informado" como erro. Isso evita viés de classificação.
    df["Prescricao_Nao_Informada"] = aud_presc.str.contains("não informado|nao informado", case=False, na=False).astype(int)
    df["Prescricao_Avaliavel"] = (~aud_presc.str.contains("não informado|nao informado|individualizada", case=False, na=False)).astype(int)
    df["Erro_Prescricao_Avaliavel"] = np.where(
        df["Prescricao_Avaliavel"].eq(1),
        aud_presc.str.contains("prolongada|curta|acima|abaixo", case=False, na=False).astype(int),
        np.nan,
    )

    # Administração: variável mais informativa no banco atual.
    df["Erro_Administracao"] = aud_admin.str.contains(
        "incompleta|suspensão precoce|suspensao precoce", case=False, na=False
    ).astype(int)
    df["Administracao_Prolongada"] = aud_admin.str.contains("prolongada", case=False, na=False).astype(int)
    df["Administracao_Completa"] = aud_admin.str.contains("ciclo completo", case=False, na=False).astype(int)
    df["Administracao_Sem_Base"] = aud_admin.str.contains("sem base", case=False, na=False).astype(int)

    df["Obito"] = desfecho.str.contains("óbito|obito", case=False, na=False).astype(int)
    df["Cultura_Sim"] = cultura.str.contains("^sim|sim", case=False, na=False).astype(int)
    df["Amplo_Espectro"] = espectro.str.contains("amplo", case=False, na=False).astype(int)
    df["Tempo_Determinado"] = tempo_det.str.contains("^sim", case=False, na=False).astype(int)

    # Garante colunas essenciais para gráficos, mesmo se vierem ausentes.
    for col in ["CUSTO_TOTAL_R$", "DOT_Exato", "Idade"]:
        if col not in df.columns:
            df[col] = np.nan
    for col in ["ATB", "CLASSE", "Sindrome", "Motivo", "Sexo", "Desfecho", "AUD_PRESCRICAO", "AUD_ADMINISTRACAO"]:
        if col not in df.columns:
            df[col] = np.nan

    return df


def agregar_internacoes(df):
    """Agrega registros de prescrição para análise por internação."""
    col_presc = nome_coluna_prescricao(df)
    presc_agg = (col_presc, "nunique") if col_presc else ("ID_Internacao", "size")

    agg = df.groupby("ID_Internacao", dropna=True).agg(
        n_prescricoes=presc_agg,
        idade=("Idade", "first"),
        sexo=("Sexo", "first"),
        setor_principal=("Setor_Padronizado", moda_segura),
        custo_total=("CUSTO_TOTAL_R$", "sum"),
        dot_total=("DOT_Exato", "sum"),
        ddd_total=("DDD_TOTAL", "sum") if "DDD_TOTAL" in df.columns else ("DOT_Exato", "sum"),
        n_atb=("ATB", "nunique"),
        amplo_espectro=("Amplo_Espectro", "max"),
        erro_admin_algum=("Erro_Administracao", "max"),
        tempo_determinado_algum=("Tempo_Determinado", "max"),
        cultura_alguma=("Cultura_Sim", "max"),
        obito=("Obito", "max"),
        desfecho_final=("Desfecho", desfecho_final),
        desfechos_registrados=("Desfecho", lambda x: " | ".join(pd.unique(x.dropna().astype(str)))),
    ).reset_index()

    agg["multiplos_desfechos"] = agg["desfechos_registrados"].astype(str).str.contains(r"\|").astype(int)
    return agg

# ============================================================
# 4. FUNÇÕES DE TESTES ESTATÍSTICOS
# ============================================================


def teste_categorico(tabela):
    """Aplica Fisher para 2x2; qui-quadrado para tabelas maiores."""
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
    """Cria tabela 2x2 garantindo linhas/colunas 0 e 1."""
    tab = pd.crosstab(df[exposicao], df[desfecho])
    tab = tab.reindex(index=[0, 1], columns=[0, 1], fill_value=0)
    return tab


def fisher_or_com_ic(tabela):
    """Fisher + OR com IC95% por statsmodels Table2x2."""
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
    """Pós-teste par-a-par com correção de Holm-Bonferroni."""
    linhas = []
    for (g1, s1), (g2, s2) in itertools.combinations(dados, 2):
        if len(s1) > 0 and len(s2) > 0:
            u, p = stats.mannwhitneyu(s1, s2, alternative="two-sided")
            linhas.append({
                "Grupo 1": g1,
                "Grupo 2": g2,
                "U de Mann-Whitney": u,
                "p bruto": p,
            })
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
    """Executa Kruskal-Wallis e pós-teste, mostrando interpretação."""
    res = kruskal_por_grupo(df, grupo, valor)
    if res is None:
        st.warning(f"Dados insuficientes para {titulo}.")
        return
    h, p, dados = res
    mostrar_resultado_teste(titulo, estatistica=h, p=p)
    if p < ALFA:
        st.write("**Pós-teste:** Mann-Whitney par-a-par com correção de Holm-Bonferroni.")
        st.dataframe(pares_mannwhitney_holm(dados), use_container_width=True)


def correlacoes_spearman(df, variaveis):
    """Gera tabela com rho, p-valor e interpretação para todos os pares."""
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

# ============================================================
# 5. LEITURA DOS DADOS
# ============================================================

uploaded = st.sidebar.file_uploader("Enviar planilha Excel ou CSV", type=["xlsx", "xls", "csv"])

try:
    df_original = carregar_arquivo(uploaded)
    df = preparar_variaveis(df_original)
    df_internacao = agregar_internacoes(df)
except Exception as e:
    st.error(f"Erro ao carregar dados: {e}")
    st.stop()

# ============================================================
# 6. FILTROS
# ============================================================

st.sidebar.header("Filtros")
setores = sorted(df["Setor_Padronizado"].dropna().unique().tolist())
setor_selecionado = st.sidebar.multiselect(
    "Setor de internação",
    options=setores,
    default=setores,
)

motivos = sorted(df["Motivo"].dropna().astype(str).unique().tolist())
motivo_selecionado = st.sidebar.multiselect(
    "Motivo da prescrição",
    options=motivos,
    default=motivos,
)

df_filtrado = df[
    df["Setor_Padronizado"].isin(setor_selecionado)
    & df["Motivo"].astype(str).isin(motivo_selecionado)
].copy()

ids_filtrados = df_filtrado["ID_Internacao"].dropna().unique().tolist()
df_internacao_filtrado = df_internacao[df_internacao["ID_Internacao"].isin(ids_filtrados)].copy()

st.caption(
    "Convenção estatística: p<0,05 sugere evidência estatística de associação/diferença; "
    "p≥0,05 indica ausência de evidência suficiente, não prova ausência de efeito. "
    "As análises inferenciais são exploratórias quando a amostra é pequena."
)

# ============================================================
# 7. ABAS DO DASHBOARD
# ============================================================

tab1, tab2, tab3, tab4, tab5, tab6, tab7, tab8 = st.tabs([
    "👁️ 1. Visão Geral",
    "📈 2. Padrões Clínicos",
    "🏥 3. Desfechos e Risco",
    "💰 4. Farmacoeconomia",
    "🎯 5. Auditoria Avançada",
    "🗄️ 6. Dados",
    "📋 7. Qualidade & Artigo",
    "🔮 8. Modelagem",
])

# ------------------------------------------------------------
# ABA 1 - VISÃO GERAL
# ------------------------------------------------------------
with tab1:
    c1, c2, c3, c4, c5 = st.columns(5)
    c1.metric("Internações", f"{df_internacao_filtrado['ID_Internacao'].nunique()}")
    c2.metric("Prescrições", f"{len(df_filtrado)}")
    c3.metric("Antimicrobianos", f"{df_filtrado['ATB'].nunique()}")
    c4.metric("Custo total", f"R$ {df_filtrado['CUSTO_TOTAL_R$'].sum():,.2f}")
    c5.metric("DOT mediano", f"{df_filtrado['DOT_Exato'].median():.2f}")

    st.markdown("---")
    col1, col2 = st.columns(2)
    with col1:
        st.plotly_chart(
            px.pie(df_filtrado, names="Setor_Padronizado", title="Prescrições por setor", hole=0.45),
            use_container_width=True,
        )
    with col2:
        top_atb = df_filtrado["ATB"].value_counts().reset_index()
        top_atb.columns = ["ATB", "n"]
        st.plotly_chart(
            px.bar(top_atb.head(15), x="n", y="ATB", orientation="h", title="Antimicrobianos mais prescritos"),
            use_container_width=True,
        )

    col3, col4 = st.columns(2)
    with col3:
        st.plotly_chart(
            px.pie(df_filtrado, names="Motivo", title="Motivo da prescrição", hole=0.45),
            use_container_width=True,
        )
    with col4:
        st.plotly_chart(
            px.pie(df_filtrado, names="Desfecho", title="Desfechos registrados nas prescrições", hole=0.45),
            use_container_width=True,
        )

# ------------------------------------------------------------
# ABA 2 - PADRÕES CLÍNICOS
# ------------------------------------------------------------
with tab2:
    st.subheader("Padrões clínicos: DDD, DOT, classes e síndromes")

    col1, col2 = st.columns(2)
    with col1:
        st.plotly_chart(
            px.pie(df_filtrado, names="Status_DDD", title="Classificação pela DDD", hole=0.45),
            use_container_width=True,
        )
    with col2:
        st.plotly_chart(
            px.box(df_filtrado, x="Setor_Padronizado", y="DOT_Exato", points="all", title="DOT por setor"),
            use_container_width=True,
        )

    testar_kruskal_e_posthoc(
        df_filtrado,
        grupo="Setor_Padronizado",
        valor="DOT_Exato",
        titulo="Kruskal-Wallis para DOT por setor",
    )

    st.markdown("---")
    col3, col4 = st.columns(2)
    with col3:
        tabela_classes = pd.crosstab(df_filtrado["Setor_Padronizado"], df_filtrado["CLASSE"])
        st.plotly_chart(
            px.bar(
                tabela_classes.reset_index(),
                x="Setor_Padronizado",
                y=tabela_classes.columns,
                barmode="stack",
                title="Classes farmacológicas por setor",
            ),
            use_container_width=True,
        )

        teste, p, efeito, menor_esp, chi2 = teste_categorico(tabela_classes)
        detalhes = f"efeito={efeito:.2f}".replace(".", ",")
        if pd.notna(menor_esp):
            detalhes += f"; menor esperado={menor_esp:.2f}".replace(".", ",")
        mostrar_resultado_teste(f"{teste} para classe farmacológica por setor", estatistica=chi2, p=p, detalhes=detalhes)
        if pd.notna(menor_esp) and menor_esp < 5:
            st.warning("Há células esperadas <5; interprete o qui-quadrado com cautela ou considere colapsar categorias.")

    with col4:
        top_sind = df_filtrado["Sindrome"].value_counts().reset_index()
        top_sind.columns = ["Síndrome", "n"]
        st.plotly_chart(
            px.bar(top_sind, x="n", y="Síndrome", orientation="h", title="Síndromes/indicações"),
            use_container_width=True,
        )

        tabela_ddd_setor = pd.crosstab(df_filtrado["Setor_Padronizado"], df_filtrado["Status_DDD"])
        st.write("Status DDD por setor")
        st.dataframe(tabela_ddd_setor, use_container_width=True)
        teste, p, efeito, menor_esp, chi2 = teste_categorico(tabela_ddd_setor)
        detalhes = f"efeito={efeito:.2f}".replace(".", ",")
        if pd.notna(menor_esp):
            detalhes += f"; menor esperado={menor_esp:.2f}".replace(".", ",")
        mostrar_resultado_teste(f"{teste} para status DDD por setor", estatistica=chi2, p=p, detalhes=detalhes)

# ------------------------------------------------------------
# ABA 3 - DESFECHOS E RISCO
# ------------------------------------------------------------
with tab3:
    st.subheader("Desfechos e risco: análise por internação")
    st.caption("Esta aba usa a base agregada por internação para reduzir pseudorrepetição.")

    col1, col2 = st.columns(2)
    with col1:
        tabela_desf = pd.crosstab(df_internacao_filtrado["setor_principal"], df_internacao_filtrado["desfecho_final"])
        st.plotly_chart(
            px.bar(
                tabela_desf.reset_index(),
                x="setor_principal",
                y=tabela_desf.columns,
                barmode="group",
                title="Desfecho final por setor principal",
            ),
            use_container_width=True,
        )

        teste, p, efeito, menor_esp, chi2 = teste_categorico(tabela_desf)
        detalhes = f"efeito={efeito:.2f}".replace(".", ",")
        if pd.notna(menor_esp):
            detalhes += f"; menor esperado={menor_esp:.2f}".replace(".", ",")
        mostrar_resultado_teste(f"{teste} para desfecho final por setor", estatistica=chi2, p=p, detalhes=detalhes)

    with col2:
        st.write("Resumo por internação")
        cols_resumo = [
            "ID_Internacao", "setor_principal", "idade", "n_prescricoes", "n_atb",
            "dot_total", "custo_total", "desfecho_final",
        ]
        st.dataframe(df_internacao_filtrado[cols_resumo], use_container_width=True)

    st.markdown("---")
    st.subheader("Odds Ratio exploratório 2x2")
    col3, col4, col5 = st.columns(3)
    with col3:
        exposicao = st.selectbox(
            "Exposição",
            ["amplo_espectro", "erro_admin_algum", "tempo_determinado_algum", "cultura_alguma"],
            index=0,
        )
    with col4:
        desfecho_bin = st.selectbox(
            "Desfecho binário",
            ["obito", "erro_admin_algum", "amplo_espectro"],
            index=0,
        )
    with col5:
        st.caption("OR bruto com Fisher exato. Não interpretar como causalidade.")

    if exposicao == desfecho_bin:
        st.warning("Escolha exposição e desfecho diferentes.")
    else:
        tab_or = tabela_2x2_binaria(df_internacao_filtrado, exposicao, desfecho_bin)
        tab_or.index = [f"{exposicao}=0", f"{exposicao}=1"]
        tab_or.columns = [f"{desfecho_bin}=0", f"{desfecho_bin}=1"]
        st.dataframe(tab_or, use_container_width=True)
        odds, li, ls, p = fisher_or_com_ic(tab_or)
        detalhes = f"OR={odds:.2f}; IC95% {li:.2f}–{ls:.2f}".replace(".", ",")
        mostrar_resultado_teste("Fisher exato para OR bruto", p=p, detalhes=detalhes)

# ------------------------------------------------------------
# ABA 4 - FARMACOECONOMIA
# ------------------------------------------------------------
with tab4:
    st.subheader("Farmacoeconomia: custo por prescrição e por internação")

    col1, col2 = st.columns(2)
    with col1:
        st.plotly_chart(
            px.box(
                df_filtrado,
                x="Setor_Padronizado",
                y="CUSTO_TOTAL_R$",
                points="all",
                title="Custo por prescrição",
            ),
            use_container_width=True,
        )
    with col2:
        st.plotly_chart(
            px.scatter(
                df_filtrado,
                x="DOT_Exato",
                y="CUSTO_TOTAL_R$",
                color="Setor_Padronizado",
                trendline="ols",
                title="Custo por DOT - prescrição",
            ),
            use_container_width=True,
        )

    testar_kruskal_e_posthoc(
        df_filtrado,
        grupo="Setor_Padronizado",
        valor="CUSTO_TOTAL_R$",
        titulo="Kruskal-Wallis para custo por setor",
    )

    st.markdown("---")
    st.subheader("Custo total por internação")
    col3, col4 = st.columns(2)
    with col3:
        st.plotly_chart(
            px.box(
                df_internacao_filtrado,
                x="setor_principal",
                y="custo_total",
                points="all",
                title="Custo total por internação",
            ),
            use_container_width=True,
        )
    with col4:
        st.plotly_chart(
            px.scatter(
                df_internacao_filtrado,
                x="dot_total",
                y="custo_total",
                color="setor_principal",
                size="n_atb",
                title="Custo total x DOT total por internação",
            ),
            use_container_width=True,
        )

    testar_kruskal_e_posthoc(
        df_internacao_filtrado,
        grupo="setor_principal",
        valor="custo_total",
        titulo="Kruskal-Wallis para custo total por internação segundo setor principal",
    )

    corr_custo = correlacoes_spearman(df_internacao_filtrado, ["dot_total", "custo_total", "n_atb", "n_prescricoes"])
    st.write("Correlação farmacoeconômica por internação")
    st.dataframe(corr_custo.round(4), use_container_width=True)

# ------------------------------------------------------------
# ABA 5 - AUDITORIA AVANÇADA
# ------------------------------------------------------------
with tab5:
    st.subheader("Auditoria avançada: associações categóricas")

    col1, col2 = st.columns(2)
    with col1:
        tb_admin_setor = pd.crosstab(df_filtrado["Setor_Padronizado"], df_filtrado["Erro_Administracao"])
        st.write("Erro/incompletude de administração por setor")
        st.dataframe(tb_admin_setor, use_container_width=True)
        teste, p, efeito, menor_esp, chi2 = teste_categorico(tb_admin_setor)
        detalhes = f"efeito={efeito:.2f}".replace(".", ",")
        if pd.notna(menor_esp):
            detalhes += f"; menor esperado={menor_esp:.2f}".replace(".", ",")
        mostrar_resultado_teste(f"{teste} para administração por setor", estatistica=chi2, p=p, detalhes=detalhes)
        if pd.notna(menor_esp) and menor_esp < 5:
            st.warning("Há células esperadas <5; interprete com cautela.")

    with col2:
        tb_admin_motivo = pd.crosstab(df_filtrado["Motivo"], df_filtrado["Erro_Administracao"])
        st.write("Erro/incompletude de administração por motivo")
        st.dataframe(tb_admin_motivo, use_container_width=True)
        teste, p, efeito, menor_esp, chi2 = teste_categorico(tb_admin_motivo)
        detalhes = f"efeito={efeito:.2f}".replace(".", ",")
        if pd.notna(menor_esp):
            detalhes += f"; menor esperado={menor_esp:.2f}".replace(".", ",")
        mostrar_resultado_teste(f"{teste} para administração por motivo", estatistica=chi2, p=p, detalhes=detalhes)
        if pd.notna(menor_esp) and menor_esp < 5:
            st.warning("Há células esperadas <5; considere colapsar categorias ou ampliar a amostra.")

    st.markdown("---")
    col3, col4 = st.columns(2)
    with col3:
        tb_ddd_admin = pd.crosstab(df_filtrado["Status_DDD"], df_filtrado["Erro_Administracao"])
        st.write("Status DDD x erro/incompletude de administração")
        st.dataframe(tb_ddd_admin, use_container_width=True)
        teste, p, efeito, menor_esp, chi2 = teste_categorico(tb_ddd_admin)
        detalhes = f"efeito={efeito:.2f}".replace(".", ",")
        if pd.notna(menor_esp):
            detalhes += f"; menor esperado={menor_esp:.2f}".replace(".", ",")
        mostrar_resultado_teste(f"{teste} para status DDD x administração", estatistica=chi2, p=p, detalhes=detalhes)

    with col4:
        tb_tempo_admin = pd.crosstab(df_filtrado["Tempo_Determinado"], df_filtrado["Erro_Administracao"])
        tb_tempo_admin.index = ["Tempo não determinado", "Tempo determinado"][:len(tb_tempo_admin.index)]
        st.write("Tempo determinado x erro/incompletude de administração")
        st.dataframe(tb_tempo_admin, use_container_width=True)
        if tb_tempo_admin.shape == (2, 2):
            odds, li, ls, p = fisher_or_com_ic(tb_tempo_admin)
            detalhes = f"OR={odds:.2f}; IC95% {li:.2f}–{ls:.2f}".replace(".", ",")
            mostrar_resultado_teste("Fisher exato para tempo determinado x administração", p=p, detalhes=detalhes)
        else:
            teste, p, efeito, menor_esp, chi2 = teste_categorico(tb_tempo_admin)
            mostrar_resultado_teste(f"{teste} para tempo determinado x administração", estatistica=chi2, p=p)

    st.markdown("---")
    st.subheader("Auditoria da prescrição: qualidade do preenchimento")
    aud = df_filtrado["AUD_PRESCRICAO"].value_counts(dropna=False).reset_index()
    aud.columns = ["Categoria", "n"]
    st.dataframe(aud, use_container_width=True)
    pct_nao_info = df_filtrado["Prescricao_Nao_Informada"].mean() * 100 if len(df_filtrado) else np.nan
    if pd.notna(pct_nao_info) and pct_nao_info > 50:
        st.warning(
            f"{pct_nao_info:.1f}% das prescrições estão como 'Não informado'. "
            "Não use AUD_PRESCRICAO como desfecho primário sem qualificar a ausência de informação."
        )

# ------------------------------------------------------------
# ABA 6 - DADOS
# ------------------------------------------------------------
with tab6:
    st.subheader("Base de prescrições filtrada")
    st.dataframe(df_filtrado, use_container_width=True, height=420)
    st.download_button(
        "Baixar prescrições tratadas em CSV",
        data=df_filtrado.to_csv(index=False, sep=";", decimal=",").encode("utf-8-sig"),
        file_name="prescricoes_tratadas.csv",
        mime="text/csv",
    )

    st.subheader("Base agregada por internação")
    st.dataframe(df_internacao_filtrado, use_container_width=True, height=320)
    st.download_button(
        "Baixar base por internação em CSV",
        data=df_internacao_filtrado.to_csv(index=False, sep=";", decimal=",").encode("utf-8-sig"),
        file_name="internacoes_agregadas.csv",
        mime="text/csv",
    )

# ------------------------------------------------------------
# ABA 7 - QUALIDADE E ARTIGO
# ------------------------------------------------------------
with tab7:
    st.subheader("Qualidade do banco e saídas para artigo")

    n_presc = len(df_filtrado)
    n_int = df_internacao_filtrado["ID_Internacao"].nunique()
    n_multi_desf = int(df_internacao_filtrado["multiplos_desfechos"].sum()) if len(df_internacao_filtrado) else 0
    n_nao_info = int(df_filtrado["Prescricao_Nao_Informada"].sum()) if n_presc else 0
    n_cultura = int(df_filtrado["Cultura_Sim"].sum()) if n_presc else 0
    n_tempo = int(df_filtrado["Tempo_Determinado"].sum()) if n_presc else 0
    n_admin_erro = int(df_filtrado["Erro_Administracao"].sum()) if n_presc else 0
    n_ddd_ok = int(df_filtrado["DDD_Adequada"].sum()) if n_presc else 0
    n_obitos = int(df_internacao_filtrado["obito"].sum()) if len(df_internacao_filtrado) else 0

    k1, k2, k3, k4, k5 = st.columns(5)
    k1.metric("Prescrições", n_presc)
    k2.metric("Internações", n_int)
    k3.metric("Óbitos", n_obitos)
    k4.metric("Múltiplos desfechos", n_multi_desf)
    k5.metric("Setores", df_filtrado["Setor_Padronizado"].nunique())

    indicadores = pd.DataFrame([
        ["Adequação à DDD", n_ddd_ok, n_presc, formatar_pct_ic(n_ddd_ok, n_presc)],
        ["Administração incompleta/suspensa", n_admin_erro, n_presc, formatar_pct_ic(n_admin_erro, n_presc)],
        ["Tempo determinado pelo prescritor", n_tempo, n_presc, formatar_pct_ic(n_tempo, n_presc)],
        ["Cultura registrada como sim", n_cultura, n_presc, formatar_pct_ic(n_cultura, n_presc)],
        ["AUD_PRESCRICAO não informado", n_nao_info, n_presc, formatar_pct_ic(n_nao_info, n_presc)],
        ["Óbito por internação", n_obitos, n_int, formatar_pct_ic(n_obitos, n_int)],
    ], columns=["Indicador", "n", "denominador", "% e IC95%"])
    st.dataframe(indicadores, use_container_width=True)

    st.markdown("---")
    st.write("Tabela 1 — descrição numérica por prescrição")
    tabela_num_presc = pd.DataFrame({
        "Idade": descrever_num(df_filtrado["Idade"]),
        "DOT_Exato": descrever_num(df_filtrado["DOT_Exato"]),
        "CUSTO_TOTAL_R$": descrever_num(df_filtrado["CUSTO_TOTAL_R$"]),
        "Ratio_DDD": descrever_num(df_filtrado["Ratio_DDD"]),
    }).T
    st.dataframe(tabela_num_presc.round(3), use_container_width=True)

    st.write("Tabela 2 — descrição numérica por internação")
    tabela_num_int = pd.DataFrame({
        "idade": descrever_num(df_internacao_filtrado["idade"]),
        "dot_total": descrever_num(df_internacao_filtrado["dot_total"]),
        "custo_total": descrever_num(df_internacao_filtrado["custo_total"]),
        "n_prescricoes": descrever_num(df_internacao_filtrado["n_prescricoes"]),
        "n_atb": descrever_num(df_internacao_filtrado["n_atb"]),
    }).T
    st.dataframe(tabela_num_int.round(3), use_container_width=True)

    st.markdown("---")
    st.write("Texto-base para método/resultados")
    st.code(
        f"""Estudo observacional retrospectivo baseado em {n_presc} prescrições antimicrobianas referentes a {n_int} internações.
As análises foram conduzidas em duas unidades: prescrição e internação. Para variáveis contínuas assimétricas, utilizaram-se mediana e intervalo interquartil, com comparação por Kruskal-Wallis ou Mann-Whitney. Proporções foram apresentadas com IC95% pelo método de Wilson. Variáveis categóricas foram avaliadas por Fisher exato em tabelas 2x2 ou qui-quadrado quando aplicável.
A adequação à DDD foi observada em {formatar_pct_ic(n_ddd_ok, n_presc)} das prescrições. Administração incompleta ou suspensão precoce ocorreu em {formatar_pct_ic(n_admin_erro, n_presc)}. Cultura registrada como sim ocorreu em {formatar_pct_ic(n_cultura, n_presc)}. A proporção de AUD_PRESCRICAO não informado foi de {formatar_pct_ic(n_nao_info, n_presc)}, indicando limitação de completude para inferência sobre adequação da prescrição médica.""",
        language="markdown",
    )

# ------------------------------------------------------------
# ABA 8 - MODELAGEM
# ------------------------------------------------------------
with tab8:
    st.subheader("Modelagem, correlação e limites inferenciais")

    st.markdown("### 1. Correlação de Spearman — prescrição")
    vars_corr_presc = ["Idade", "DOT_Exato", "CUSTO_TOTAL_R$", "Ratio_DDD"]
    df_corr_presc = df_filtrado[vars_corr_presc].dropna()

    col1, col2 = st.columns([1, 2])
    if len(df_corr_presc) > 5:
        rho, pmat = stats.spearmanr(df_corr_presc)
        rho_df = pd.DataFrame(rho, index=vars_corr_presc, columns=vars_corr_presc)
        p_df = pd.DataFrame(pmat, index=vars_corr_presc, columns=vars_corr_presc)
        with col2:
            fig = px.imshow(
                rho_df,
                text_auto=".2f",
                color_continuous_scale="RdBu_r",
                zmin=-1,
                zmax=1,
                title="Matriz de correlação de Spearman",
            )
            st.plotly_chart(fig, use_container_width=True)
        with col1:
            st.write("Matriz de p-valores")
            st.dataframe(p_df.round(4), use_container_width=True)
            st.caption("Correlação não implica causalidade; use como análise exploratória.")

        interp_corr = correlacoes_spearman(df_filtrado, vars_corr_presc)
        st.write("Interpretação par-a-par")
        st.dataframe(interp_corr.round(4), use_container_width=True)
    else:
        st.warning("Dados insuficientes para correlação de Spearman por prescrição.")

    st.markdown("---")
    st.markdown("### 2. Correlação de Spearman — internação")
    vars_corr_int = ["idade", "dot_total", "custo_total", "n_prescricoes", "n_atb"]
    interp_corr_int = correlacoes_spearman(df_internacao_filtrado, vars_corr_int)
    st.dataframe(interp_corr_int.round(4), use_container_width=True)

    st.markdown("---")
    st.markdown("### 3. Mortalidade — análise exploratória")
    eventos = int(df_internacao_filtrado["obito"].sum()) if len(df_internacao_filtrado) else 0
    n_modelo = len(df_internacao_filtrado)
    st.write(f"Eventos de óbito por internação: **{eventos}/{n_modelo}**.")

    if eventos < 10:
        st.warning(
            "Não é metodologicamente recomendado ajustar regressão logística múltipla com menos de 10 eventos. "
            "A saída abaixo deve ser descrita como exploratória/descritiva."
        )
        comparacoes = []
        for var in ["idade", "dot_total", "custo_total", "n_prescricoes", "n_atb"]:
            g0 = pd.to_numeric(df_internacao_filtrado.loc[df_internacao_filtrado["obito"] == 0, var], errors="coerce").dropna()
            g1 = pd.to_numeric(df_internacao_filtrado.loc[df_internacao_filtrado["obito"] == 1, var], errors="coerce").dropna()
            if len(g0) > 0 and len(g1) > 0:
                u, p = stats.mannwhitneyu(g1, g0, alternative="two-sided")
                comparacoes.append({
                    "Variável": var,
                    "Mediana óbito": g1.median(),
                    "Mediana não óbito": g0.median(),
                    "U de Mann-Whitney": u,
                    "p-valor": p,
                    "p formatado": formatar_p(p),
                    "Decisão": classe_decisao_p(p),
                    "Interpretação": interpretar_p(p),
                })
        st.dataframe(pd.DataFrame(comparacoes).round(4), use_container_width=True)
    else:
        base = df_internacao_filtrado[["obito", "idade", "amplo_espectro", "erro_admin_algum", "n_atb"]].dropna()
        X = sm.add_constant(base[["idade", "amplo_espectro", "erro_admin_algum", "n_atb"]])
        y = base["obito"]
        try:
            modelo = sm.Logit(y, X).fit(disp=False)
            res = pd.DataFrame({
                "OR ajustado": np.exp(modelo.params),
                "IC95% inferior": np.exp(modelo.conf_int()[0]),
                "IC95% superior": np.exp(modelo.conf_int()[1]),
                "p-valor": modelo.pvalues,
            })
            res["p formatado"] = res["p-valor"].apply(formatar_p)
            res["Decisão"] = res["p-valor"].apply(classe_decisao_p)
            res["Interpretação"] = res["p-valor"].apply(interpretar_p)
            st.dataframe(res.drop("const", errors="ignore"), use_container_width=True)
        except Exception as e:
            st.warning(f"O modelo não convergiu: {e}")

    st.markdown("---")
    st.markdown("### 4. Análise de sobrevivência")
    st.warning(
        "Kaplan-Meier não deve ser interpretado como sobrevida hospitalar nesta base enquanto não houver "
        "data real de admissão, alta, óbito ou tempo até evento. O DOT mede exposição terapêutica, não tempo de sobrevida."
    )

    st.markdown("---")
    st.markdown("### 5. Conclusão estatística automática")
    st.info(
        "Para a redação do artigo, priorize DOT, DDD, custos, incompletude da administração e qualidade do registro. "
        "Regressões sobre óbito devem ser evitadas enquanto houver poucos eventos. Use p-values sempre acompanhados de tamanho de efeito, IC95% quando aplicável e ressalva de análise exploratória."
    )
