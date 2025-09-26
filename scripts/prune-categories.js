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
]);

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
  let removedCount = 0;
  let changedWords = 0;
  const result = [];

  for (const entry of entries) {
    if (entry && Array.isArray(entry.categories)) {
      const before = entry.categories.length;
      const filtered = entry.categories.filter((c) => !CATEGORIES_TO_REMOVE.has(String(c)));
      removedCount += before - filtered.length;
      if (filtered.length !== before) changedWords++;

      if (dropEmpty && filtered.length === 0) {
        continue; // remove this entry entirely
      }
      result.push({ ...entry, categories: filtered });
    } else {
      // Not an expected shape; pass through untouched
      result.push(entry);
    }
  }

  return { data: result, removedCount, changedWords };
}

function pruneObjectModel(map, dropEmpty) {
  let removedCount = 0;
  let changedWords = 0;
  const out = {};

  for (const [word, categories] of Object.entries(map)) {
    if (Array.isArray(categories)) {
      const before = categories.length;
      const filtered = categories.filter((c) => !CATEGORIES_TO_REMOVE.has(String(c)));
      removedCount += before - filtered.length;
      if (filtered.length !== before) changedWords++;
      if (dropEmpty && filtered.length === 0) {
        continue; // drop key entirely
      }
      out[word] = filtered;
    } else {
      out[word] = categories;
    }
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


