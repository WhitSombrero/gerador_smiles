"""Geração de SMILES base: anéis decorados com auxocromos aleatórios."""

import random


def _inserir_fragmento_em_posicao(smiles: str, posicao: int, fragmento: str) -> str:
    """Insere `(fragmento)` na posição especificada do SMILES."""
    posicao_valida = max(0, min(posicao, len(smiles)))
    return smiles[:posicao_valida] + f"({fragmento})" + smiles[posicao_valida:]


def _adicionar_auxocromos(
    smiles_base: str,
    auxocromos: dict[int, str],
    posicoes_disponiveis: list[int],
    quantidade: int,
) -> str:
    """
    Insere `quantidade` auxocromos aleatórios em `smiles_base`.

    As inserções são feitas da maior para a menor posição
    para não deslocar os índices anteriores.
    """
    posicoes_escolhidas = random.sample(posicoes_disponiveis, quantidade)
    insercoes = [
        (pos, auxocromos[random.randint(1, len(auxocromos))])
        for pos in posicoes_escolhidas
    ]
    insercoes.sort(reverse=True)  # maior posição primeiro

    smiles = smiles_base
    for posicao, frag in insercoes:
        smiles = _inserir_fragmento_em_posicao(smiles, posicao, frag)
    return smiles


def gerar_fragmentos_com_auxocromos(
    aneis: dict[int, str],
    auxocromos: dict[int, str],
    repeticoes: int = 200,
    posicoes_disponiveis: list[int] | None = None,
    quantidade_auxocromos: int = 3,
) -> dict[int, str]:
    """
    Gera variantes de cada anel com auxocromos inseridos aleatoriamente.

    Parâmetros
    ----------
    aneis : dict[int, str]
        Dicionário de SMILES base dos anéis.
    auxocromos : dict[int, str]
        Dicionário de fragmentos substituintes disponíveis.
    repeticoes : int
        Quantidade de variantes por anel.
    posicoes_disponiveis : list[int] | None
        Posições do SMILES onde auxocromos podem ser inseridos.
        Padrão: [2, 3, 4, 5, 6].
    quantidade_auxocromos : int
        Quantos auxocromos inserir por variante.

    Retorna
    -------
    dict[int, str]
        Dicionário {id_sequencial: smiles_gerado}.
    """
    if posicoes_disponiveis is None:
        posicoes_disponiveis = [2, 3, 4, 5, 6]

    resultado: dict[int, str] = {}
    contador = 1

    for smiles_base in aneis.values():
        for _ in range(repeticoes):
            smiles_gerado = _adicionar_auxocromos(
                smiles_base, auxocromos, posicoes_disponiveis, quantidade_auxocromos
            )
            resultado[contador] = smiles_gerado
            contador += 1

    return resultado
