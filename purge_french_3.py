import os
import re

strings_to_translate = {
    "Variantes modifiées : Le modèle COMPASS précédent a été purgé.": "Modified variants: Previous COMPASS model purged.",
    "2/2 - MCMC en arrière-plan. La session est débloquée, vous pouvez lancer d'autres analyses.": "2/2 - Background MCMC. Session unlocked, you can run other analyses.",
    "COMPASS a filtré toutes les cellules (doublets purs).": "COMPASS filtered all cells (pure doublets).",
    "COMPASS terminé : Architecture clonale et Matrice Multi-Omique verrouillées.": "COMPASS finished: Clonal architecture and Multi-Omics Matrix locked.",
    "Calcul UMAP et évaluations terminés.": "UMAP computation and evaluations finished."
}

french_comment_patterns = [
    "pulvérise", "extrémités", "ignore", "logiques",
    "contrôleur d'affichage", "déclenche", "bloqué", "contrôles & axes", "remplacement du",
    "quand une heatmap", "quand tu renommes", "quand compass recalcule", "contourne les bugs", "lié",
    "nettoyage strict", "vire le bruit technique", "déroulante", "coloration",
    "isolement", "tout le reste devient", "coloré par-dessus", "respecte aussi la convention",
    "moteur de rendu récursif", "plongée récursive", "registre local à la session",
    "pointe vers l'env", "vérifier qu'il y a quelque chose à sauver", "alignement strict",
    "effacer proprement", "bloc de calcul réactif", "bleu transparent élégant",
    "interface épurée", "calcul du top marqueurs", "défensif", "destructeur a été atomisé",
    "proprement à partir", "actualise le menu", "dessiné en premier", "le plus bas",
    "jointure sécurisée", "échantillonnage stochastique", "topologie globale",
    "exactes de", "vectorisé", "masques booléens", "remplace la boucle",
    "logique de révélation", "à adapter", "selon la structure", "exemple : si le type", "hypothétique basé",
    "agrégation statistique", "discret et élégant", "remise à zéro",
    "peut cibler un plot précis", "prendront leurs valeurs", "selon ton thème",
    "calcul immédiat", "clones purs", "récupère systématiquement", "vérité terrain",
    "rafraîchissement", "détruit pas le filtre", "tournera ici jusqu'à",
    "renvoie l'image", "récepteur silencieux", "injecte les cibles", "inférence data-driven",
    "indépendante de", "plus de na", "calcul est terminé", "contient l'arbre",
    "ne s'affiche que si", "gestionnaire de téléchargement", "crashs de librairies",
    "on met à jour", "on veut que", "on s'assure", "on génère", "on y injecte",
    "on ne met à jour", "on ne détruit", "on exclut", "on contourne"
]

def process_file(filepath):
    with open(filepath, 'r', encoding='utf-8') as f:
        lines = f.readlines()
        
    new_lines = []
    for line in lines:
        # String translations
        for fr, en in strings_to_translate.items():
            if fr in line:
                line = line.replace(fr, en)
        
        # Detect remaining French comments
        is_french_comment = False
        if '#' in line:
            comment_part = line[line.find('#'):].lower()
            # Also check for accents in comments
            if re.search(r'[éèêëàâäîïôöùûüç]', comment_part):
                is_french_comment = True
            for w in french_comment_patterns:
                if w in comment_part:
                    is_french_comment = True
                    break
                    
        if is_french_comment:
            if line[:line.find('#')].strip() == '':
                continue
            else:
                line = line[:line.find('#')].rstrip() + '\n'
                
        new_lines.append(line)
        
    with open(filepath, 'w', encoding='utf-8') as f:
        f.writelines(new_lines)

for root, dirs, files in os.walk('R'):
    for f in files:
        if f.endswith('.R'):
            process_file(os.path.join(root, f))

print("Phase 3 purge completed.")
