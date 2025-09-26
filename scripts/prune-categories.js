'use strict';

const fs = require('fs');
const path = require('path');

const DEFAULT_INPUT = path.join(__dirname, '..', 'data', 'words.json');
const DEFAULT_OUTPUT = path.join(__dirname, '..', 'data', 'words.pruned.json');

const CATEGORIES_TO_REMOVE = new Set([
  '21st Century',
  '20th Century',
  '2020s',
  '2010s',
  'Things American',
  "Flower-class Corvettes"
]);

// Also remove categories that begin with these phrases
const CATEGORY_PREFIX_PATTERNS = [
  /^Starts with\b/,
  /^Ends with\b/,
];

function shouldRemoveCategory(category) {
  const name = String(category);
  if (CATEGORIES_TO_REMOVE.has(name)) return true;
  for (const re of CATEGORY_PREFIX_PATTERNS) {
    if (re.test(name)) return true;
  }
  return false;
}

function parseArgs(argv) {
  const args = Object.create(null);
  for (let i = 2; i < argv.length; i++) {
    const token = argv[i];
    if (token.startsWith('--')) {
      const [key, maybeVal] = token.split('=');
      const k = key.replace(/^--/, '');
      if (maybeVal !== undefined) {
        args[k] = maybeVal;
      } else {
        const next = argv[i + 1];
        if (next && !next.startsWith('--')) {
          args[k] = next; i++;
        } else {
          args[k] = true;
        }
      }
    }
  }
  return args;
}

function pruneArrayModel(entries, dropEmpty) {
  // Normalize and apply explicit/prefix removals first
  const items = [];
  for (const entry of entries) {
    if (entry && Array.isArray(entry.categories)) {
      const originalLen = entry.categories.length;
      const filtered = entry.categories.filter((c) => !shouldRemoveCategory(c));
      items.push({ entry, categories: filtered, originalLen });
    } else {
      items.push({ entry, categories: undefined, originalLen: 0 });
    }
  }

  function recomputeCounts(arr) {
    const counts = new Map();
    for (const it of arr) {
      if (!Array.isArray(it.categories)) continue;
      for (const c of it.categories) counts.set(c, (counts.get(c) || 0) + 1);
    }
    return counts;
  }

  // Iteratively remove categories that appear only once across all entries
  while (true) {
    const counts = recomputeCounts(items);
    const singles = new Set(Array.from(counts.entries()).filter(([, n]) => n === 1).map(([c]) => c));
    if (singles.size === 0) break;
    for (const it of items) {
      if (!Array.isArray(it.categories)) continue;
      it.categories = it.categories.filter((c) => !singles.has(c));
    }
  }

  // Build final list and compute stats
  let removedCount = 0;
  let changedWords = 0;
  const result = [];
  for (const it of items) {
    if (!Array.isArray(it.categories)) {
      result.push(it.entry);
      continue;
    }
    const finalLen = it.categories.length;
    removedCount += Math.max(0, it.originalLen - finalLen);
    if (it.originalLen !== finalLen) changedWords++;
    if (dropEmpty && finalLen === 0) {
      continue;
    }
    result.push({ ...it.entry, categories: it.categories });
  }

  return { data: result, removedCount, changedWords };
}

function pruneObjectModel(map, dropEmpty) {
  // Apply explicit/prefix-based removals first
  const working = {};
  const originalLens = {};
  for (const [word, categories] of Object.entries(map)) {
    if (Array.isArray(categories)) {
      originalLens[word] = categories.length;
      working[word] = categories.filter((c) => !shouldRemoveCategory(c));
    } else {
      originalLens[word] = 0;
      working[word] = categories;
    }
  }

  function recomputeCounts(obj) {
    const counts = new Map();
    for (const cats of Object.values(obj)) {
      if (!Array.isArray(cats)) continue;
      for (const c of cats) counts.set(c, (counts.get(c) || 0) + 1);
    }
    return counts;
  }

  // Iteratively remove categories that appear only once across all words
  let current = working;
  while (true) {
    const counts = recomputeCounts(current);
    const singles = new Set(Array.from(counts.entries()).filter(([, n]) => n === 1).map(([c]) => c));
    if (singles.size === 0) break;
    const next = {};
    for (const [w, cats] of Object.entries(current)) {
      if (!Array.isArray(cats)) { next[w] = cats; continue; }
      const newCats = cats.filter((c) => !singles.has(c));
      if (dropEmpty && newCats.length === 0) {
        continue; // drop word entirely
      }
      next[w] = newCats;
    }
    current = next;
  }

  const out = current;
  let removedCount = 0;
  let changedWords = 0;
  for (const [w, origLen] of Object.entries(originalLens)) {
    const finalCats = out.hasOwnProperty(w) && Array.isArray(out[w]) ? out[w] : [];
    const finalLen = finalCats.length;
    removedCount += Math.max(0, origLen - finalLen);
    if (origLen !== finalLen || !out.hasOwnProperty(w)) changedWords++;
  }

  return { data: out, removedCount, changedWords };
}

function detectModel(json) {
  if (Array.isArray(json)) return 'array';
  if (json && typeof json === 'object') return 'object';
  return 'unknown';
}

function main() {
  const args = parseArgs(process.argv);
  const inPlace = Boolean(args['in-place'] || args.inplace);
  const dropEmpty = Boolean(args['drop-empty'] || args.dropEmpty);
  const inputPath = args.input ? path.resolve(process.cwd(), args.input) : DEFAULT_INPUT;
  const outputPath = inPlace
    ? inputPath
    : (args.output ? path.resolve(process.cwd(), args.output) : DEFAULT_OUTPUT);

  if (!fs.existsSync(inputPath)) {
    console.error(`Input file not found: ${inputPath}`);
    process.exit(1);
  }

  const raw = fs.readFileSync(inputPath, 'utf8');
  let json;
  try {
    json = JSON.parse(raw);
  } catch (err) {
    console.error('Failed to parse JSON from', inputPath);
    console.error(err.message);
    process.exit(1);
  }

  const model = detectModel(json);
  let result;
  if (model === 'array') {
    result = pruneArrayModel(json, dropEmpty);
  } else if (model === 'object') {
    result = pruneObjectModel(json, dropEmpty);
  } else {
    console.error('Unsupported JSON structure. Expected array or object at root.');
    process.exit(1);
  }

  const { data, removedCount, changedWords } = result;

  // If writing in-place and the source looked minified (no newlines), keep minified
  const looksMinified = !/\n/.test(raw);
  const spacing = looksMinified && inPlace ? undefined : 2;
  const serialized = JSON.stringify(data, null, spacing) + (spacing ? '\n' : '');

  // If in-place, create a backup first
  if (inPlace) {
    const backupPath = inputPath + '.bak';
    try {
      fs.copyFileSync(inputPath, backupPath, fs.constants.COPYFILE_EXCL);
      console.log(`Backup written: ${backupPath}`);
    } catch (e) {
      // If backup exists, continue without failing
      if (e && e.code !== 'EEXIST') throw e;
    }
  }

  fs.writeFileSync(outputPath, serialized, 'utf8');

  console.log(`Removed ${removedCount} category tags across ${changedWords} words.`);
  console.log(`Output written to: ${outputPath}`);
  if (!inPlace) {
    console.log('Tip: use --in-place to overwrite the input file, with a .bak backup.');
  }
}

if (require.main === module) {
  try {
    main();
  } catch (err) {
    console.error(err);
    process.exit(1);
  }
}


