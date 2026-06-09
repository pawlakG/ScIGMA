import os

def process_file(filepath):
    with open(filepath, 'r', encoding='utf-8') as f:
        lines = f.readlines()
        
    new_lines = []
    for line in lines:
        if '# UPDATED' not in line.upper():
            new_lines.append(line)
            
    with open(filepath, 'w', encoding='utf-8') as f:
        f.writelines(new_lines)

for root, dirs, files in os.walk('R'):
    for f in files:
        if f.endswith('.R'):
            process_file(os.path.join(root, f))

print("Purge completed.")
