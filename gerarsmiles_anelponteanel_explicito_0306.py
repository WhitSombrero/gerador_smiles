
# ============================================================
# 2. Função de Fusão e Rastreador de Colunas
# ============================================================
def gerar_complexos(df1: pd.DataFrame, df2: pd.DataFrame, pontes_dict: dict, n_combinacoes: int) -> pd.DataFrame:
    dados_complexos = []
    chaves_pontes_disponiveis = list(pontes_dict.keys())

    for _ in range(n_combinacoes):
        row1 = df1.sample(1).iloc[0]
        row2 = df2.sample(1).iloc[0]

        chave_ponte = random.choice(chaves_pontes_disponiveis)
        ponte_str = pontes_dict[chave_ponte]

        smi_anel1 = row1['SMILES_Modificado']
        smi_anel2 = row2['SMILES_Modificado']
        smi_complexo = f"{smi_anel1}{ponte_str}{smi_anel2}"

        linha = {}

        for col in df1.columns:
            linha[f"{col}_Anel1"] = row1[col]

        linha["ID_Ponte"] = str(chave_ponte)
        linha["Ponte_Estrutura"] = ponte_str

        for col in df2.columns:
            linha[f"{col}_Anel2"] = row2[col]

        linha["SMILES_Complexo"] = smi_complexo
        dados_complexos.append(linha)

    return pd.DataFrame(dados_complexos)

# ============================================================
# 3. Limpeza RDKit
# ============================================================
def limpar_dataframe_complexos(df: pd.DataFrame) -> pd.DataFrame:
    total_inicial = len(df)

    df_limpo = df.drop_duplicates(subset=['SMILES_Complexo']).copy()
    total_apos_duplicatas = len(df_limpo)

    def complexo_eh_valido(smi):
        mol = Chem.MolFromSmiles(smi)
        return mol is not None

    mascara_validos = df_limpo['SMILES_Complexo'].apply(complexo_eh_valido)
    df_limpo = df_limpo[mascara_validos].copy()
    total_apos_rdkit = len(df_limpo)

    print("\n--- RELATÓRIO DE FUSÃO E LIMPEZA ---")
    print(f"Fusões tentadas:           {total_inicial}")
    print(f"Removidas por duplicata:   {total_inicial - total_apos_duplicatas}")
    print(f"Removidas por valência:    {total_apos_duplicatas - total_apos_rdkit}")
    print(f"TOTAL COMPLEXOS VÁLIDOS:   {total_apos_rdkit}")

    return df_limpo.reset_index(drop=True)

# ============================================================
# 4. Interface Interativa (Seleção Simplificada)
# ============================================================
def iniciar_fusao():
    print("\n=============================================")
    print("      CONSTRUTOR DE COMPLEXOS (FUSÃO)      ")
    print("=============================================\n")

    # 1. Escaneia a memória por DataFrames
    variaveis_globais = globals()
    dfs_disponiveis = {}
    contador = 1

    for nome, var in variaveis_globais.items():
        # Filtra variáveis que são DataFrames e começam com "df_"
        if isinstance(var, pd.DataFrame) and nome.startswith("df_"):
            dfs_disponiveis[contador] = nome
            contador += 1

    if not dfs_disponiveis:
        print("[ERRO] Nenhum DataFrame (Anel) foi encontrado na memória.")
        print("Certifique-se de rodar a geração de anéis primeiro!")
        return None, None

    # 2. Mostra o menu simplificado
    print("Anéis disponíveis para fusão:")
    for num, nome in dfs_disponiveis.items():
        print(f"  [{num}] {nome}")

    try:
        escolha1 = int(input("\nDigite o NÚMERO do DataFrame para o Anel 1: "))
        escolha2 = int(input("Digite o NÚMERO do DataFrame para o Anel 2: "))
    except ValueError:
        print("\n[ERRO] Por favor, digite apenas o número.")
        return None, None

    if escolha1 not in dfs_disponiveis or escolha2 not in dfs_disponiveis:
        print("\n[ERRO] Escolha inválida. Número não encontrado na lista.")
        return None, None

    df1_carregado = variaveis_globais[dfs_disponiveis[escolha1]]
    df2_carregado = variaveis_globais[dfs_disponiveis[escolha2]]

    # 3. Seleção de Pontes e Quantidade
    print(f"\nBancos de Pontes disponíveis: {list(CATALOGO_PONTES_GLOBAL.keys())}")
    nome_pontes = input("Digite o nome do banco de PONTES: ").strip()

    if nome_pontes not in CATALOGO_PONTES_GLOBAL:
        print("\n[ERRO] Dicionário de pontes não encontrado.")
        return None, None

    try:
        qtd_fusoes = int(input("\nQuantas fusões (combinações aleatórias) gerar? "))
    except ValueError:
        print("\n[ERRO] Por favor, digite um número inteiro válido.")
        return None, None

    # 4. Geração e Limpeza
    print("\nRealizando cruzamentos e conectando as pontes... aguarde.")
    pontes_escolhidas = CATALOGO_PONTES_GLOBAL[nome_pontes]
    df_complexos_bruto = gerar_complexos(df1_carregado, df2_carregado, pontes_escolhidas, qtd_fusoes)
    df_complexos_limpo = limpar_dataframe_complexos(df_complexos_bruto)

    # 5. Nomeando o Output Simplificado
    nome_padrao = f"df_complexo_{escolha1}x{escolha2}"
    print(f"\nNome sugerido para o DataFrame final: {nome_padrao}")
    nome_personalizado = input("Digite um nome curto (ou pressione Enter para aceitar a sugestão): ").strip()

    if nome_personalizado:
        nome_df_complexo = f"df_{nome_personalizado}" if not nome_personalizado.startswith("df_") else nome_personalizado
    else:
        nome_df_complexo = nome_padrao

    return df_complexos_limpo, nome_df_complexo

# ============================================================
# PONTO DE ENTRADA DA FUSÃO
# ============================================================
if __name__ == "__main__":
    resultado_fusao = iniciar_fusao()

    if resultado_fusao[0] is not None:
        df_complexo_final, nome_dinamico_complexo = resultado_fusao

        if not df_complexo_final.empty:
            globals()[nome_dinamico_complexo] = df_complexo_final

            print(f"\n✅ Variável salva com sucesso: '{nome_dinamico_complexo}'")
            print("\nVisualização das primeiras linhas:")
            pd.set_option('display.max_columns', None)
            print(globals()[nome_dinamico_complexo].head())