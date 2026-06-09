import os
import re
import json

french_indicators = [
    r'\bpour\b', r'\bavec\b', r'\bdans\b', r'\bfichier\b', r'\butilisateurs?\b',
    r'\berreur\b', r'\bajout\b', r'\bmise\b', r'\bjour\b', r'\bsi\b', 
    r'\bça\b', r'\boui\b', r'\bnon\b', r'\bici\b', r'\bbien\b', r'\btrès\b',
    r'\bfait\b', r'\btout\b', r'\brien\b', r'\bcet\b', r'\bcette\b', r'\bces\b',
    r'\bceci\b', r'\bcela\b', r'\bqui\b', r'\bque\b', r'\bquoi\b', r'\bdont\b',
    r'\boù\b', r'\bquand\b', r'\bcomment\b', r'\bpourquoi\b', r'\bou\b',
    r'[éèêëàâäîïôöùûüç]' # Any French accent
]

pattern = re.compile('|'.join(french_indicators), re.IGNORECASE)

results = {}

for root, dirs, files in os.walk('R'):
    for f in files:
        if f.endswith('.R'):
            path = os.path.join(root, f)
            with open(path, 'r', encoding='utf-8') as file:
                lines = file.readlines()
            
            for i, line in enumerate(lines):
                # Ignore lines starting with roxygen2 comments or pure code without strings/comments
                if line.strip().startswith("#'"):
                    continue
                
                # Check for comments or strings
                is_match = False
                
                if '#' in line:
                    comment_part = line[line.find('#'):]
                    if pattern.search(comment_part):
                        is_match = True
                        
                # Check for strings (crude regex for strings)
                strings = re.findall(r'"([^"]*)"|\'([^\']*)\'', line)
                for s_tuple in strings:
                    s = s_tuple[0] if s_tuple[0] else s_tuple[1]
                    if pattern.search(s):
                        # Filter out common english matches like "Pour", "Si", etc.
                        is_match = True
                        
                if is_match:
                    if path not in results:
                        results[path] = []
                    results[path].append({"line": i+1, "content": line.strip()})

print(json.dumps(results, indent=2, ensure_ascii=False))
