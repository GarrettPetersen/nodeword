'use strict';

const fs = require('fs');
const path = require('path');

const WORDS_PATH = path.join(__dirname, '..', 'data', 'words.json');
const EMOJIS_PATH = path.join(__dirname, '..', 'data', 'category_emojis.json');

const DEFAULT_EMOJI = '✅';

function readJson(p) {
  const raw = fs.readFileSync(p, 'utf8');
  return JSON.parse(raw);
}

function getAllCategoriesFromWords(words) {
  const set = new Set();
  if (Array.isArray(words)) {
    for (const entry of words) {
      if (entry && Array.isArray(entry.categories)) {
        for (const c of entry.categories) set.add(String(c));
      }
    }
    return set;
  }
  // object model { word: [categories] }
  for (const cats of Object.values(words)) {
    if (Array.isArray(cats)) for (const c of cats) set.add(String(c));
  }
  return set;
}

function main() {
  if (!fs.existsSync(WORDS_PATH)) {
    console.error('words.json not found at', WORDS_PATH);
    process.exit(1);
  }
  if (!fs.existsSync(EMOJIS_PATH)) {
    console.error('category_emojis.json not found at', EMOJIS_PATH);
    process.exit(1);
  }

  const words = readJson(WORDS_PATH);
  const categoriesInUse = getAllCategoriesFromWords(words);
  const emojis = readJson(EMOJIS_PATH);

  const out = {};
  let removed = 0;
  let added = 0;

  // Add existing mappings only if category is used
  for (const [cat, emoji] of Object.entries(emojis)) {
    if (categoriesInUse.has(cat)) {
      out[cat] = emoji;
    } else {
      removed++;
    }
  }

  // Add missing categories with default emoji
  for (const cat of categoriesInUse) {
    if (!out.hasOwnProperty(cat)) {
      out[cat] = DEFAULT_EMOJI;
      added++;
    }
  }

  // Preserve pretty formatting
  const serialized = JSON.stringify(out, null, 2) + '\n';
  fs.writeFileSync(EMOJIS_PATH, serialized, 'utf8');

  console.log(`Synced categories: removed ${removed} unused; added ${added} missing.`);
  console.log(`Updated file: ${EMOJIS_PATH}`);
}

if (require.main === module) {
  try {
    main();
  } catch (err) {
    console.error(err);
    process.exit(1);
  }
}


