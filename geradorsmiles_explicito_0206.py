"""
Módulo responsável pela geração de SMILES modificados a partir de catálogos de anéis e fragmentos,
incluindo validação estrutural via RDKit.
"""

import random
import pandas as pd
from rdkit import Chem

# ============================================================
# 1. Funções de Processamento
# ============================================================

def inserir_fragmento(smiles: str, pos: int, frag: str) -> str:
    """Insere um fragmento em uma posição específica de uma string SMILES."""
    pos = max(0, min(pos, len(smiles)))
    return smiles[:pos] + f"({frag})" + smiles[pos:]

def gerar_smiles_modificado(smiles_base: str, frags_dict: dict, posicoes_possiveis: list, n_fragmentos: int):
    """Sorteia posições e fragmentos para gerar uma nova string SMILES."""
    smiles_mod = smiles_base
    n_fragmentos = min(n_fragmentos, len(posicoes_possiveis))

    posicoes = random.sample(posicoes_possiveis, n_fragmentos)
    chaves_frags_disponiveis = list(frags_dict.keys())
    chaves_sorteadas = random.choices(chaves_frags_disponiveis, k=n_fragmentos)
    fragmentos = [frags_dict[k] for k in chaves_sorteadas]

    insercoes = list(zip(posicoes, fragmentos, chaves_sorteadas))
    insercoes.sort(key=lambda x: x[0], reverse=True)

    posicoes_usadas = []
    frags_usados_str = []
    chaves_frags_str = []

    for pos, frag, chave in insercoes:
        smiles_mod = inserir_fragmento(smiles_mod, pos, frag)
        posicoes_usadas.append(pos)
        frags_usados_str.append(frag)
        chaves_frags_str.append(str(chave))

    posicoes_usadas.reverse()
    frags_usados_str.reverse()
    chaves_frags_str.reverse()

    return smiles_mod, posicoes_usadas, frags_usados_str, chaves_frags_str

def gerar_dataframe_smiles(aneis_dict: dict, frags_dict: dict, n_fragmentos: int = 1, repeticoes: int = 10) -> pd.DataFrame:
    """Gera um DataFrame contendo todas as combinações de SMILES criadas."""
    dados = []
    posicoes_possiveis = [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16]

    for id_anel, smiles_base in aneis_dict.items():
        for _ in range(repeticoes):
            smiles_mod, posicoes_usadas, frags_usados, chaves_frags = gerar_smiles_modificado(
                smiles_base,
                frags_dict,
                posicoes_possiveis,
                n_fragmentos
            )

            linha_dados = {
                "ID_Anel": id_anel,
                "SMILES_Base": smiles_base,
            }

            for i in range(len(posicoes_usadas)):
                linha_dados[f"ID_Frag_{i+1}"] = chaves_frags[i]
                linha_dados[f"Frag_Estrutura_{i+1}"] = frags_usados[i]
                linha_dados[f"Posicao_{i+1}"] = posicoes_usadas[i]

            linha_dados["SMILES_Modificado"] = smiles_mod
            dados.append(linha_dados)

    return pd.DataFrame(dados)

# ============================================================
# 2. Função de Limpeza (RDKit Validation + Duplicates)
# ============================================================

def limpar_dataframe_gerado(df: pd.DataFrame) -> pd.DataFrame:
    """Remove duplicatas e filtra apenas moléculas quimicamente válidas via RDKit."""
    total_inicial = len(df)

    df_limpo = df.drop_duplicates(subset=['SMILES_Modificado']).copy()
    total_apos_duplicatas = len(df_limpo)

    def smiles_eh_valido(smi):
        mol = Chem.MolFromSmiles(smi)
        return mol is not None

    mascara_validos = df_limpo['SMILES_Modificado'].apply(smiles_eh_valido)
    df_limpo = df_limpo[mascara_validos].copy()
    total_apos_rdkit = len(df_limpo)

    print("\n--- RELATÓRIO DE LIMPEZA DE DADOS ---")
    print(f"Total gerado inicialmente: {total_inicial}")
    print(f"Removidas por duplicata:   {total_inicial - total_apos_duplicatas}")
    print(f"Removidas por SMILES inv.: {total_apos_duplicatas - total_apos_rdkit}")
    print(f"Total final válido:        {total_apos_rdkit}\n")

    return df_limpo.reset_index(drop=True)

# ============================================================
# 3. Pipeline Principal
# ============================================================

def pipeline_gerar_e_limpar(aneis_dict: dict, frags_dict: dict, n_fragmentos: int, repeticoes: int) -> pd.DataFrame:
    """
    Executa o fluxo completo: cruza os dicionários, gera as moléculas e limpa os dados.
    Retorna o DataFrame final validado.
    """
    print("Processando cruzamento e inserções... aguarde.")
    df_bruto = gerar_dataframe_smiles(
        aneis_dict=aneis_dict,
        frags_dict=frags_dict,
        n_fragmentos=n_fragmentos,
        repeticoes=repeticoes
    )

    print("Limpando e validando moléculas com RDKit... isso pode levar alguns segundos.")
    df_resultado_limpo = limpar_dataframe_gerado(df_bruto)
    
    return df_resultado_limpo