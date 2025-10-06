// Seeded RNG utilities for deterministic daily puzzles
let CURRENT_RNG = { random: () => Math.random() };
function xmur3(str) { // hash to 32-bit seed
  let h = 1779033703 ^ str.length;
  for (let i = 0; i < str.length; i++) {
    h = Math.imul(h ^ str.charCodeAt(i), 3432918353);
    h = (h << 13) | (h >>> 19);
  }
  return function() {
    h = Math.imul(h ^ (h >>> 16), 2246822507);
    h = Math.imul(h ^ (h >>> 13), 3266489909);
    return (h ^= h >>> 16) >>> 0;
  };
}
function mulberry32(a) {
  return function() {
    let t = (a += 0x6D2B79F5);
    t = Math.imul(t ^ (t >>> 15), t | 1);
    t ^= t + Math.imul(t ^ (t >>> 7), t | 61);
    return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
  };
}
function makeRng(seedString) {
  const seed = xmur3(String(seedString))();
  const rand = mulberry32(seed);
  return {
    random: rand,
    int(max) { return Math.floor(rand() * max); },
    choice(arr) { return arr[this.int(arr.length)]; },
    shuffle(arr) {
      for (let i = arr.length - 1; i > 0; i--) {
        const j = Math.floor(rand() * (i + 1));
        [arr[i], arr[j]] = [arr[j], arr[i]];
      }
      return arr;
    }
  };
}
function withSeed(seedString, fn) {
  const prev = CURRENT_RNG;
  let result;
  CURRENT_RNG = makeRng(seedString);
  try {
    result = fn();
    if (result && typeof result.then === 'function') {
      // Async: restore RNG after the promise settles
      return result.finally(() => { CURRENT_RNG = prev; });
    }
    return result;
  } finally {
    // Sync: restore immediately; for async, finally runs now but we guard above
    if (!result || typeof result.then !== 'function') CURRENT_RNG = prev;
  }
}

function getRandomInt(max) {
  return Math.floor(CURRENT_RNG.random() * max);
}

function choice(array) {
  if (!array || array.length === 0) return undefined;
  return array[getRandomInt(array.length)];
}

function shuffle(array) {
  for (let i = array.length - 1; i > 0; i--) {
    const j = Math.floor(CURRENT_RNG.random() * (i + 1));
    [array[i], array[j]] = [array[j], array[i]];
  }
  return array;
}

function buildIndexes(wordToCategories) {
  const categoryToWords = new Map();
  for (const [word, categories] of Object.entries(wordToCategories)) {
    for (const category of categories) {
      const list = categoryToWords.get(category) || new Set();
      list.add(word);
      categoryToWords.set(category, list);
    }
  }
  return categoryToWords;
}

function generatePuzzleGraph(wordToCategories, targetWordCount = 12, maxDegree = 4, maxAttempts = 5000) {
  const categoryToWords = buildIndexes(wordToCategories);
  const allWords = Object.keys(wordToCategories);
  if (allWords.length === 0) throw new Error("No words available");

  // Graph structures
  const wordSet = new Set();
  const categorySet = new Set();
  const wordDegree = new Map();
  const categoryDegree = new Map();
  const edges = new Set(); // key: `${word}||${category}`
  // Adjacency maps for current puzzle
  const currentWordToCats = new Map(); // word -> Set(categories)
  const currentCatToWords = new Map(); // category -> Set(words)

  function removeEdge(word, category) {
    const key = word + "||" + category;
    if (!edges.has(key)) return;
    edges.delete(key);
    wordDegree.set(word, (wordDegree.get(word) || 1) - 1);
    categoryDegree.set(category, (categoryDegree.get(category) || 1) - 1);
    const wc = currentWordToCats.get(word);
    if (wc) {
      wc.delete(category);
      if (wc.size === 0) currentWordToCats.delete(word);
    }
    const cw = currentCatToWords.get(category);
    if (cw) {
      cw.delete(word);
      if (cw.size === 0) currentCatToWords.delete(category);
    }
    if (!currentCatToWords.has(category)) categorySet.delete(category);
    if (!currentWordToCats.has(word)) wordSet.delete(word);
  }

  function addEdge(word, category) {
    const key = word + "||" + category;
    if (edges.has(key)) return true;
    if ((wordDegree.get(word) || 0) >= maxDegree) return false;
    if ((categoryDegree.get(category) || 0) >= maxDegree) return false;

    // New rule: prevent any pair of words from sharing more than one category in this puzzle
    const peers = currentCatToWords.get(category);
    if (peers && peers.size > 0) {
      const wordCats = currentWordToCats.get(word);
      if (wordCats && wordCats.size > 0) {
        for (const peer of peers) {
          const peerCats = currentWordToCats.get(peer);
          if (!peerCats || peerCats.size === 0) continue;
          // If word and peer already share any category, adding this category would make it 2
          const [small, large] = wordCats.size < peerCats.size ? [wordCats, peerCats] : [peerCats, wordCats];
          for (const c of small) {
            if (large.has(c)) {
              return false;
            }
          }
        }
      }
    }

    edges.add(key);
    wordSet.add(word);
    categorySet.add(category);
    wordDegree.set(word, (wordDegree.get(word) || 0) + 1);
    categoryDegree.set(category, (categoryDegree.get(category) || 0) + 1);
    // Update adjacency
    if (!currentWordToCats.has(word)) currentWordToCats.set(word, new Set());
    currentWordToCats.get(word).add(category);
    if (!currentCatToWords.has(category)) currentCatToWords.set(category, new Set());
    currentCatToWords.get(category).add(word);
    return true;
  }

  function getWordNeighbors(word) {
    const cats = wordToCategories[word] || [];
    return cats.filter((c) => edges.has(word + "||" + c));
  }

  function getCategoryNeighbors(category) {
    const words = Array.from(categoryToWords.get(category) || []);
    return words.filter((w) => edges.has(w + "||" + category));
  }

  function isContiguous() {
    if (wordSet.size === 0) return true;
    // BFS across bipartite graph
    const visitedW = new Set();
    const visitedC = new Set();
    const startWord = wordSet.values().next().value;
    const queue = [{ type: 'word', id: startWord }];
    visitedW.add(startWord);
    while (queue.length) {
      const node = queue.shift();
      if (node.type === 'word') {
        for (const category of (wordToCategories[node.id] || [])) {
          if (edges.has(node.id + "||" + category) && !visitedC.has(category)) {
            visitedC.add(category);
            queue.push({ type: 'category', id: category });
          }
        }
      } else {
        const candidates = Array.from(categoryToWords.get(node.id) || []);
        for (const word of candidates) {
          if (edges.has(word + "||" + node.id) && !visitedW.has(word)) {
            visitedW.add(word);
            queue.push({ type: 'word', id: word });
          }
        }
      }
    }
    return visitedW.size === wordSet.size; // we only require words to be connected via categories
  }

  function violatesPairRule() {
    // Ensure no two words share more than one category within the current puzzle
    const words = Array.from(wordSet);
    for (let i = 0; i < words.length; i++) {
      const wi = words[i];
      const ci = currentWordToCats.get(wi) || new Set();
      if (ci.size <= 1) continue;
      for (let j = i + 1; j < words.length; j++) {
        const wj = words[j];
        const cj = currentWordToCats.get(wj) || new Set();
        if (cj.size <= 1) continue;
        // Count intersection quickly
        let shared = 0;
        const [small, large] = ci.size < cj.size ? [ci, cj] : [cj, ci];
        for (const c of small) {
          if (large.has(c)) {
            shared++;
            if (shared > 1) return true;
          }
        }
      }
    }
    return false;
  }

  function violatesCoverageRule() {
    // For each category present in the puzzle, ensure that every selected word
    // that has this category is actually connected to it in the graph.
    for (const category of categorySet) {
      const connectedWords = new Set(getCategoryNeighbors(category));
      for (const word of wordSet) {
        const cats = wordToCategories[word] || [];
        if (cats.includes(category)) {
          if (!connectedWords.has(word)) {
            return true; // word should connect to category but doesn't
          }
        }
      }
    }
    return false;
  }

  function attemptAddWordWithRequiredCategories(newWord, requiredCatsArray) {
    const newlyAdded = [];
    const requiredList = Array.from(new Set(requiredCatsArray || []));
    const wordCats = new Set(wordToCategories[newWord] || []);
    // Build the full set of categories this word must connect to: required + any existing puzzle categories it belongs to
    const needed = new Set(requiredList);
    for (const c of categorySet) {
      if (wordCats.has(c)) needed.add(c);
    }
    // Pre-check capacity on the word
    const wordCapacity = maxDegree - (wordDegree.get(newWord) || 0);
    if (needed.size > wordCapacity) return false;
    // Pre-check capacity on each category
    for (const c of needed) {
      if ((categoryDegree.get(c) || 0) >= maxDegree) return false;
    }
    // Add edges in order: required first, then the rest
    const ordered = requiredList.concat(Array.from(needed).filter(c => !requiredList.includes(c)));
    for (const c of ordered) {
      if (!addEdge(newWord, c)) {
        // rollback
        for (let i = newlyAdded.length - 1; i >= 0; i--) {
          removeEdge(newWord, newlyAdded[i]);
        }
        return false;
      }
      newlyAdded.push(c);
    }
    return newlyAdded;
  }

  function ensureCategoryCoverage(category, excludeWord) {
    const added = [];
    const connectees = [];
    for (const w of wordSet) {
      if (w === excludeWord) continue;
      const cats = wordToCategories[w] || [];
      if (cats.includes(category) && !edges.has(w + "||" + category)) {
        connectees.push(w);
      }
    }
    for (const w of connectees) {
      if (!addEdge(w, category)) {
        // rollback
        for (let i = added.length - 1; i >= 0; i--) {
          removeEdge(added[i], category);
        }
        return false;
      }
      added.push(w);
    }
    return true;
  }

  for (let attempt = 0; attempt < maxAttempts; attempt++) {
    // Reset graph state each attempt
    wordSet.clear(); categorySet.clear(); edges.clear(); wordDegree.clear(); categoryDegree.clear();
    currentWordToCats.clear(); currentCatToWords.clear();

    const startWord = choice(allWords);
    const startCats = (wordToCategories[startWord] || []).slice();
    if (!startCats || startCats.length === 0) continue;
    shuffle(startCats);
    let seeded = false;
    for (const firstCategory of startCats) {
      if (!addEdge(startWord, firstCategory)) { continue; }
      const candidates = Array.from(categoryToWords.get(firstCategory) || []);
      shuffle(candidates);
      let addedSecond = false;
      for (const w of candidates) {
        if (w === startWord) continue;
        if (addEdge(w, firstCategory)) { addedSecond = true; break; }
      }
      if (addedSecond) { seeded = true; break; }
      // rollback and try another category for the start word
      removeEdge(startWord, firstCategory);
    }
    if (!seeded) continue;

    // Greedily grow the graph while keeping degrees <= maxDegree and maintaining connectivity
    while (wordSet.size < targetWordCount) {
      // Special strategy for the final word: try to add a word that connects to 2 categories
      if (wordSet.size === targetWordCount - 1) {
        const candidateWords = shuffle(allWords.filter(w => !wordSet.has(w)).slice());
        let addedFinal = false;
        for (const w of candidateWords) {
          const cats = (wordToCategories[w] || []).slice();
          if (cats.length < 2) continue;
          // Prefer categories already in the graph to improve connectivity
          const catScore = (c) => (categorySet.has(c) ? 1 : 0);
          cats.sort((a, b) => catScore(b) - catScore(a));
          // Try pairs of categories; ensure at least one is already present
          for (let i = 0; i < cats.length && !addedFinal; i++) {
            for (let j = i + 1; j < cats.length; j++) {
              const c1 = cats[i], c2 = cats[j];
              if (!categoryToWords.has(c1) || !categoryToWords.has(c2)) continue;
              if (!categorySet.has(c1) && !categorySet.has(c2)) continue; // keep connectivity
          // Try to add the word with both categories (and any other applicable puzzle categories)
          if (attemptAddWordWithRequiredCategories(w, [c1, c2])) {
            addedFinal = true;
          }
              break;
            }
          }
          if (addedFinal) break;
        }
        if (addedFinal) {
          // Completed graph with enhanced connectivity
          break;
        }
        // If none found, fall back to normal growth below
      }

      const pickExistingWord = (CURRENT_RNG.random() < 0.5) || wordSet.size === 1;
      const anchorWord = pickExistingWord ? choice(Array.from(wordSet)) : choice(allWords);
      const availableCategories = (wordToCategories[anchorWord] || []).slice();
      if (availableCategories.length === 0) break;
      shuffle(availableCategories);

      let expanded = false;
      for (const category of availableCategories) {
        // ensure category exists in data
        const candidateWords = Array.from(categoryToWords.get(category) || []);
        if (candidateWords.length === 0) continue;
        shuffle(candidateWords);

        // Try to connect anchor -> category, then category -> newWord
        if (!addEdge(anchorWord, category)) continue;

        let foundNewWord = false;
        for (const w of candidateWords) {
          if (wordSet.has(w)) {
            // Also add the edge if within degree cap to increase density
            if (addEdge(w, category)) {
              foundNewWord = true;
              // we added a redundant within-graph edge, keep searching for a new word
            }
            continue;
          }
          if (attemptAddWordWithRequiredCategories(w, [category])) {
            foundNewWord = true;
            expanded = true;
            break;
          }
        }

        if (foundNewWord) break;
        // rollback edge anchor->category if it didn't lead to expansion
        removeEdge(anchorWord, category);
      }

      if (!expanded) break;
      if (!isContiguous()) break;
    }

    if (wordSet.size >= targetWordCount && isContiguous() && !violatesPairRule() && !violatesCoverageRule()) {
      // Build final structure
      const words = Array.from(wordSet);
      const categories = Array.from(categorySet);
      const edgesList = [];
      for (const key of edges) {
        const [w, c] = key.split('||');
        edgesList.push({ word: w, category: c });
      }

      // Post-step: prune orphan categories (degree 1) after placing last word
      const catCounts = new Map();
      for (const e of edgesList) {
        catCounts.set(e.category, (catCounts.get(e.category) || 0) + 1);
      }
      const prunedEdges = edgesList.filter(e => (catCounts.get(e.category) || 0) > 1);
      const prunedCategories = Array.from(new Set(prunedEdges.map(e => e.category)));
      // Keep the same words list; they should still have other connections due to earlier constraints
      return { words, categories: prunedCategories, edges: prunedEdges };
    }
  }

  throw new Error("Failed to generate a contiguous graph within constraints");
}

async function fetchWordData() {
  console.log('[Nodeword] Fetching words.json…');
  let res = await fetch('/data/words.json');
  if (!res.ok) {
    console.warn('[Nodeword] Fetch /data/words.json failed with', res.status, '— retrying relative path');
    res = await fetch('data/words.json');
  }
  if (!res.ok) throw new Error('Failed to load words.json');
  const data = await res.json();
  const totalWords = Object.keys(data || {}).length;
  console.log('[Nodeword] Loaded words.json with', totalWords, 'entries');
  return data;
}

async function fetchCategoryEmojis() {
  try {
    let res = await fetch('/data/category_emojis.json');
    if (!res.ok) res = await fetch('data/category_emojis.json');
    if (!res.ok) throw new Error('Failed to load category_emojis.json');
    const data = await res.json();
    return data;
  } catch (e) {
    console.warn('[Nodeword] category_emojis.json not loaded:', e.message);
    return {};
  }
}

const NODEWORD_CONFIG = {
  targetWords: 12,
  maxDegree: 4,
  maxAttempts: 6000,
};

function toAliasGraph(graph) {
  const wordAliases = new Map();
  const catAliases = new Map();
  const nodes = [];

  graph.words.forEach((w, i) => {
    const alias = `W${i + 1}`;
    wordAliases.set(w, alias);
    nodes.push({ id: w, alias, type: 'word' });
  });
  graph.categories.forEach((c, i) => {
    const alias = `C${i + 1}`;
    catAliases.set(c, alias);
    nodes.push({ id: c, alias, type: 'category' });
  });

  const links = graph.edges.map(e => ({
    source: wordAliases.get(e.word),
    target: catAliases.get(e.category)
  }));

  return { nodes, links };
}

function measureWordCircle(word) {
  // Fit based on the longest token; wrap other tokens if helpful but size by max token
  const text = String(word).trim();
  const tokens = text.split(/\s+/);
  const longest = tokens.reduce((m, t) => Math.max(m, t.length), 0);
  const minR = 18;
  const maxR = 72;
  const minFont = 9;
  const maxFont = 18;

  // Radius driven primarily by the longest word
  let r = Math.min(maxR, Math.max(minR, 14 + Math.ceil(longest * 1.8)));
  let font = Math.min(maxFont, Math.max(minFont, Math.floor(r * 0.44)));

  function layout(radius, fontPx) {
    const maxWidthPx = radius * 1.7;
    const avgCharPx = fontPx * 0.58;
    const maxChars = Math.max(3, Math.floor(maxWidthPx / avgCharPx));
    const maxLines = Math.max(1, Math.floor((radius * 1.6) / (fontPx + 2)));
    const lines = [];
    let line = '';
    for (const tk of tokens) {
      const next = line ? line + ' ' + tk : tk;
      if (next.length <= maxChars) line = next; else { if (line) lines.push(line); line = tk; }
    }
    if (line) lines.push(line);
    while (lines.length > maxLines) {
      let idx = 0; let len = 0;
      for (let i = 0; i < lines.length; i++) { if (lines[i].length > len) { len = lines[i].length; idx = i; } }
      const ws = lines[idx].split(' ');
      if (ws.length <= 1) break;
      const half = Math.ceil(ws.length / 2);
      lines.splice(idx, 1, ws.slice(0, half).join(' '), ws.slice(half).join(' '));
    }
    const totalH = lines.length * (fontPx + 2);
    const fits = (totalH <= radius * 1.6) && (Math.max(...lines.map(l => l.length), 0) * avgCharPx <= maxWidthPx);
    return { lines, totalH, fits };
  }

  for (let i = 0; i < 24; i++) {
    const attempt = layout(r, font);
    if (attempt.fits) {
      const startY = -attempt.totalH / 2 + font / 2;
      return { radius: r, fontSize: font, parts: attempt.lines, startY };
    }
    // Prefer increasing radius first (up to cap), then shrinking font
    if (r < maxR) r += 4; else if (font > minFont) font -= 1; else break;
  }
  const startY = -font / 2;
  return { radius: r, fontSize: font, parts: [text], startY };
}

function renderForceGraph(container, aliasGraph, wordToCategories, categoryEmojis, persist) {
  container.innerHTML = '';
  const svg = d3.select(container).append('svg').attr('class', 'graph-svg');
  let width = container.clientWidth;
  let height = Math.min(820, Math.max(520, Math.floor(window.innerHeight * 0.75)));
  svg.attr('width', width).attr('height', height);
  const scene = svg.append('g').attr('class', 'scene');

  const nodes = aliasGraph.nodes.map(n => ({ ...n }));
  const links = aliasGraph.links.map(l => ({ ...l }));
  const aliasToNode = new Map(nodes.map(n => [n.alias, n]));

  console.log('[Nodeword] Rendering force graph:', {
    words: nodes.filter(n => n.type === 'word').length,
    categories: nodes.filter(n => n.type === 'category').length,
    edges: links.length
  });

  // Prelayout with Kamada-Kawai (more iterations, early stop, random restarts) to reduce crossings
  (function kkPrelayout() {
    if (persist && persist.restore && persist.restore.nodePositions) {
      const pos = persist.restore.nodePositions;
      const map = new Map(Object.entries(pos));
      for (const n of nodes) {
        const p = map.get(n.alias);
        if (p && typeof p.x === 'number' && typeof p.y === 'number') {
          n.x = p.x; n.y = p.y; n.vx = 0; n.vy = 0;
        }
      }
      return;
    }
    const N = nodes.length;
    if (N === 0) return;
    const aliasToIndex = new Map(nodes.map((n, i) => [n.alias, i]));
    // Build adjacency for unweighted shortest paths
    const adj = Array.from({ length: N }, () => []);
    for (const e of links) {
      const si = aliasToIndex.get(e.source);
      const ti = aliasToIndex.get(e.target);
      if (si == null || ti == null) continue;
      adj[si].push(ti);
      adj[ti].push(si);
    }
    // All-pairs shortest path distances via BFS
    const INF = 1e9;
    const dist = Array.from({ length: N }, () => Array(N).fill(INF));
    for (let s = 0; s < N; s++) {
      const q = [s];
      dist[s][s] = 0;
      for (let qi = 0; qi < q.length; qi++) {
        const u = q[qi];
        for (const v of adj[u]) if (dist[s][v] === INF) {
          dist[s][v] = dist[s][u] + 1;
          q.push(v);
        }
      }
    }
    // Determine desired lengths
    let dmax = 0;
    for (let i = 0; i < N; i++) for (let j = i + 1; j < N; j++) if (dist[i][j] < INF) dmax = Math.max(dmax, dist[i][j]);
    if (dmax === 0) dmax = 1;
    const pad = Math.max(24, Math.round(Math.min(width, height) * 0.06));
    const Lmin = 40, Lmax = Math.min(width, height) / 2 - pad;
    const L = (i, j) => (dist[i][j] >= INF ? Lmax : (Lmin + (Lmax - Lmin) * (dist[i][j] / dmax)));
    const K = (i, j) => (dist[i][j] <= 0 || dist[i][j] >= INF ? 0 : 1 / (dist[i][j] * dist[i][j]));
    let bestStress = Infinity; let bestX = null, bestY = null;
    const restarts = 7;
    const kScale = 4e-5; // scale spring constants to stabilize further
    for (let rs = 0; rs < restarts; rs++) {
      // Init positions on a circle, alternating words and categories around the ring
      let x = new Array(N), y = new Array(N);
      const radius0 = Math.min(width, height) * 0.25;
      const wordsIdx = nodes.map((n, i) => ({ n, i })).filter(v => v.n.type === 'word').map(v => v.i);
      const catsIdx = nodes.map((n, i) => ({ n, i })).filter(v => v.n.type === 'category').map(v => v.i);
      const orderIdx = [];
      const maxLen = Math.max(wordsIdx.length, catsIdx.length);
      for (let k = 0; k < maxLen; k++) {
        if (k < wordsIdx.length) orderIdx.push(wordsIdx[k]);
        if (k < catsIdx.length) orderIdx.push(catsIdx[k]);
      }
      // If counts are very uneven, ensure we still place all nodes
      for (let i = 0; i < N; i++) if (!orderIdx.includes(i)) orderIdx.push(i);
      for (let p = 0; p < N; p++) {
        const idx = orderIdx[p];
        const a = (2 * Math.PI * p) / N;
        x[idx] = width / 2 + Math.cos(a) * radius0;
        y[idx] = height / 2 + Math.sin(a) * radius0;
      }
      // Gradient descent on stress
      const iters = 2500;
      const step0 = 0.0015;
      let calmIters = 0;
      for (let it = 0; it < iters; it++) {
        const fx = new Array(N).fill(0);
        const fy = new Array(N).fill(0);
        for (let i = 0; i < N; i++) {
          for (let j = i + 1; j < N; j++) {
            const dx = x[i] - x[j];
            const dy = y[i] - y[j];
            const r = Math.hypot(dx, dy) || 1e-6;
            const l = L(i, j);
            const k = K(i, j) * kScale;
            if (k === 0) continue;
            const f = k * (r - l) / r;
            const fxij = f * dx;
            const fyij = f * dy;
            fx[i] -= fxij; fy[i] -= fyij;
            fx[j] += fxij; fy[j] += fyij;
          }
        }
        let maxDelta = 0;
        const step = Math.max(0.0003, step0 * (1 - it / iters));
        for (let i = 0; i < N; i++) {
          const dx = step * fx[i];
          const dy = step * fy[i];
          x[i] -= dx; y[i] -= dy;
          const md = Math.max(Math.abs(dx), Math.abs(dy));
          if (md > maxDelta) maxDelta = md;
        }
        if (maxDelta < 0.15) {
          calmIters++;
          if (calmIters > 12) break;
        } else {
          calmIters = 0;
        }
      }
      // Compute stress
      let stress = 0;
      for (let i = 0; i < N; i++) {
        for (let j = i + 1; j < N; j++) {
          const dx = x[i] - x[j];
          const dy = y[i] - y[j];
          const r = Math.hypot(dx, dy) || 1e-6;
          const l = L(i, j);
          const k = K(i, j);
          if (k === 0) continue;
          const d = r - l; stress += k * d * d;
        }
      }
      if (stress < bestStress) { bestStress = stress; bestX = x; bestY = y; }
    }
    const x = bestX, y = bestY;
    console.log('[Nodeword] KK prelayout complete. stress=', Number(bestStress.toFixed(2)));
    // Normalize to padded viewport
    let minX = Math.min(...x), maxX = Math.max(...x);
    let minY = Math.min(...y), maxY = Math.max(...y);
    let cw = Math.max(1, maxX - minX), ch = Math.max(1, maxY - minY);
    // Target lengths range
    const Lmin2 = 80, Lmax2 = Math.min(width, height) * 0.45;
    const sx = (width - 2 * pad) / cw, sy = (height - 2 * pad) / ch;
    const s = Math.min(sx, sy);
    for (let i = 0; i < N; i++) {
      const nx = (x[i] - minX - cw / 2) * s + width / 2;
      const ny = (y[i] - minY - ch / 2) * s + height / 2;
      nodes[i].x = nx; nodes[i].y = ny; nodes[i].vx = 0; nodes[i].vy = 0;
    }
  })();

  // Set of categories that exist in THIS puzzle (restrict emoji logic to these)
  const puzzleCategories = new Set(nodes.filter(n => n.type === 'category').map(n => n.id));

  const link = scene.append('g')
    .attr('stroke-linecap', 'round')
    .selectAll('line')
    .data(links)
    .enter()
    .append('line')
    .attr('class', 'link');

  const node = scene.append('g')
    .selectAll('g')
    .data(nodes, d => d.alias)
    .enter()
    .append('g')
    .attr('class', d => `node ${d.type}`);

  // Prepare category selection (base nodes); words have an overlay for dragging
  const catNodes = node.filter(d => d.type === 'category');

  // Make categories draggable via base nodes; words are dragged via overlay to avoid double-drag
  const dragBehavior = d3.drag()
    .on('start', (event, d) => {
      if (!event.active) simulation.alphaTarget(BACKGROUND_ALPHA).restart();
      d.fx = d.x; d.fy = d.y;
    })
    .on('drag', (event, d) => {
      d.fx = event.x; d.fy = event.y;
    })
    .on('end', (event, d) => {
      d.fx = null; d.fy = null;
      if (!event.active) simulation.alphaTarget(BACKGROUND_ALPHA);
    });

  // Drag only category base nodes; words are dragged via overlay to avoid double-drag
  catNodes.call(dragBehavior);
  // Save assignment on diamond drag end as well
  catNodes.on('end.save', function() { try { snapshotAssignmentOnly(); } catch {} });

  const radius = d => d.type === 'category' ? 16 : 12;

  // Draw shapes: word nodes as circles, category nodes as diamonds (rotated squares)
  node.filter(d => d.type === 'word')
    .append('circle')
    .attr('r', radius);
  catNodes
    .append('rect')
    .attr('x', -12)
    .attr('y', -12)
    .attr('width', 24)
    .attr('height', 24)
    .attr('transform', 'rotate(45)');

  // Optional alias labels only for words
  node.filter(d => d.type === 'word')
    .append('text')
    .text(d => d.alias);

  // Emoji overlay for categories
  const catEmoji = catNodes
    .append('text')
    .attr('class', 'cat-emoji')
    .attr('text-anchor', 'middle')
    .attr('dominant-baseline', 'central')
    .text('');
  // Hidden label to show on solve
  const catLabel = catNodes
    .append('text')
    .attr('class', 'cat-label')
    .attr('text-anchor', 'middle')
    .attr('y', -20)
    .text('');

  // Build word-circle overlay assignments with minimization of pre-solved categories
  const wordNodes = nodes.filter(n => n.type === 'word');
  const wordLabelsBase = aliasGraph.nodes.filter(n => n.type === 'word').map(n => n.id);

  // Precompute connected word-aliases per category-alias using alias links
  const catToWordAliases = new Map();
  for (const l of links) {
    const catAlias = l.target; // category alias (string at this point)
    const wordAlias = l.source; // word alias
    const arr = catToWordAliases.get(catAlias) || [];
    arr.push(wordAlias);
    catToWordAliases.set(catAlias, arr);
  }

  // Custom force: pull category diamonds toward the centroid of their attached words
  function categoryCentroidForce(strength = 0.04) {
    return () => {
      for (const [catAlias, wordAliases] of catToWordAliases.entries()) {
        if (!wordAliases || wordAliases.length === 0) continue;
        const cn = aliasToNode.get(catAlias);
        if (!cn) continue;
        let cx = 0, cy = 0, n = 0;
        for (const wa of wordAliases) {
          const wn = aliasToNode.get(wa);
          if (!wn || typeof wn.x !== 'number' || typeof wn.y !== 'number') continue;
          cx += wn.x; cy += wn.y; n++;
        }
        if (n === 0) continue;
        cx /= n; cy /= n;
        cn.vx += (cx - cn.x) * strength;
        cn.vy += (cy - cn.y) * strength;
      }
    };
  }

  // Custom force: repel intersecting edges (anti-crossing)
  function antiCrossingForce(strength = 0.12) {
    function toNode(end) {
      if (!end) return null;
      if (typeof end === 'string') return aliasToNode.get(end) || null;
      // d3 may set source/target to node objects
      if (end.alias) return end;
      return null;
    }
    function segsIntersect(a1, a2, b1, b2) {
      const d = (p, q, r) => (q.x - p.x) * (r.y - p.y) - (q.y - p.y) * (r.x - p.x);
      const o1 = d(a1, a2, b1);
      const o2 = d(a1, a2, b2);
      const o3 = d(b1, b2, a1);
      const o4 = d(b1, b2, a2);
      return (o1 * o2 < 0) && (o3 * o4 < 0);
    }
    return () => {
      const L = links.length;
      for (let i = 0; i < L; i++) {
        const l1 = links[i];
        const s1 = toNode(l1.source); const t1 = toNode(l1.target);
        if (!s1 || !t1) continue;
        for (let j = i + 1; j < L; j++) {
          const l2 = links[j];
          // Skip if sharing endpoints
          if (l1.source === l2.source || l1.source === l2.target || l1.target === l2.source || l1.target === l2.target) continue;
          const s2 = toNode(l2.source); const t2 = toNode(l2.target);
          if (!s2 || !t2) continue;
          if (!segsIntersect(s1, t1, s2, t2)) continue;
          // Compute midpoints and apply perpendicular impulses to separate
          const m1x = (s1.x + t1.x) / 2, m1y = (s1.y + t1.y) / 2;
          const m2x = (s2.x + t2.x) / 2, m2y = (s2.y + t2.y) / 2;
          let vx = m1x - m2x, vy = m1y - m2y;
          let len = Math.hypot(vx, vy);
          if (len < 1e-6) { // pick a perpendicular to l1
            vx = t1.y - s1.y; vy = -(t1.x - s1.x); len = Math.hypot(vx, vy) || 1; 
          }
          vx /= len; vy /= len;
          const k = strength;
          s1.vx += vx * k; s1.vy += vy * k; t1.vx += vx * k; t1.vy += vy * k;
          s2.vx -= vx * k; s2.vy -= vy * k; t2.vx -= vx * k; t2.vy -= vy * k;
        }
      }
    };
  }

  function countSolvedForAssignment(assignMap) {
    let solved = 0;
    for (const [catAlias, wordAliases] of catToWordAliases.entries()) {
      if (!wordAliases || wordAliases.length === 0) continue;
      let intersection = null;
      for (const wa of wordAliases) {
        const w = assignMap.get(wa);
        const catsArr = wordToCategories[w] || [];
        const cats = new Set(catsArr);
        if (!intersection) {
          intersection = cats;
        } else {
          const next = new Set();
          for (const c of intersection) if (cats.has(c)) next.add(c);
          intersection = next;
        }
        if (!intersection || intersection.size === 0) break;
      }
      if (!intersection || intersection.size === 0) continue;
      // Restrict to categories in this puzzle
      let hasPuzzleCat = false;
      for (const c of intersection) { if (puzzleCategories.has(c)) { hasPuzzleCat = true; break; } }
      if (hasPuzzleCat) solved++;
    }
    return solved;
  }

  let assignment;
  const restoredSolved = Boolean(persist && persist.restore && persist.restore.solved);
  if (restoredSolved) {
    // Force canonical assignment on solved restore to present a clean solved state
    assignment = new Map();
    const orderedWordNodes = aliasGraph.nodes
      .filter(n => n.type === 'word')
      .sort((a, b) => (parseInt(String(a.alias).slice(1)) || 0) - (parseInt(String(b.alias).slice(1)) || 0));
    for (const wn of orderedWordNodes) {
      assignment.set(wn.alias, wn.id);
    }
  } else if (persist && persist.restore && persist.restore.assignment) {
    assignment = new Map(Object.entries(persist.restore.assignment));
  } else {
    let bestAssignment = null;
    let bestSolved = Infinity;
    for (let attempt = 0; attempt < 10; attempt++) {
      const labels = shuffle(wordLabelsBase.slice());
      const assign = new Map();
      for (let i = 0; i < wordNodes.length; i++) {
        assign.set(wordNodes[i].alias, labels[i % labels.length]);
      }
      const solvedCount = countSolvedForAssignment(assign);
      if (solvedCount < bestSolved) {
        bestSolved = solvedCount;
        bestAssignment = assign;
        if (bestSolved === 0) break;
      }
    }
    assignment = bestAssignment || new Map();
  }

  // Precompute circle metrics and per-node extra repulsion based on circle size
  const circleMetrics = new Map(); // node.alias -> {radius,fontSize,parts}
  for (const wn of wordNodes) {
    const w = assignment.get(wn.alias);
    const m = measureWordCircle(w);
    circleMetrics.set(wn.alias, m);
  }

  const circlesLayer = scene.append('g').attr('class', 'word-circles');
  const wordCircleG = circlesLayer.selectAll('g')
    .data(wordNodes)
    .enter()
    .append('g')
    .attr('class', 'word-circle');

  wordCircleG.append('circle')
    .attr('r', d => (circleMetrics.get(d.alias)?.radius || 24));

  // Apply immediate highlight before any interactions
  try {
    const isFinalLevelInit = Boolean(persist && persist.meta && persist.meta.levelToday === 5);
    if (isFinalLevelInit) {
      const firstWordInit = aliasGraph.nodes.filter(n=>n.type==='word')[0]?.id;
      if (firstWordInit) {
        let targetAliasInit = null;
        // assignment may not yet match nodeAliasToWord; build a temp map from initial assignment
        const initialMap = new Map(assignment);
        initialMap.forEach((w, alias) => { if (w === firstWordInit) targetAliasInit = alias; });
        if (targetAliasInit) {
          wordCircleG.classed('highlighted', false);
          wordCircleG.filter(d => d.alias === targetAliasInit).classed('highlighted', true);
        }
      }
    }
  } catch {}

  wordCircleG.each(function(d) {
    const g = d3.select(this);
    const w = assignment.get(d.alias);
    const m = circleMetrics.get(d.alias);
    for (let i = 0; i < m.parts.length; i++) {
      g.append('text')
        .attr('y', m.startY + i * (m.fontSize + 2))
        .style('font-size', m.fontSize + 'px')
        .text(m.parts[i]);
    }
  });

  // Highlight the word circle that currently holds the first puzzle word on final daily level (5)
  try {
    const isFinalLevel = Boolean(persist && persist.meta && persist.meta.levelToday === 5);
    if (isFinalLevel) {
      const firstWord = aliasGraph.nodes.filter(n=>n.type==='word')[0]?.id;
      if (firstWord) {
        let targetAlias = null;
        nodeAliasToWord.forEach((w, alias) => { if (w === firstWord) targetAlias = alias; });
        if (targetAlias) {
          wordCircleG.classed('highlighted', false);
          wordCircleG.filter(d => d.alias === targetAlias).classed('highlighted', true);
        }
      }
    }
  } catch {}

  // Click handling for swapping circles between nodes
  let selected = null; // d (node datum) for selected word-circle
  let solved = restoredSolved;

  function deselect() {
    selected = null;
    wordCircleG.classed('selected', false);
  }

  function flyTo(gSel, x, y) {
    gSel.transition().duration(280).ease(d3.easeCubicOut)
      .attr('transform', `translate(${x},${y})`);
  }

  // Make word circles draggable (overlay), with a small movement threshold so taps/clicks reliably select
  const DRAG_THRESHOLD = 8; // pixels
  const dragWordOverlay = d3.drag().clickDistance(DRAG_THRESHOLD)
    .on('start', function(event, d) {
      d._dragStartX = event.x;
      d._dragStartY = event.y;
      d._didDrag = false;
      // If nothing is selected yet, dragging can select this circle; otherwise keep current selection for swap
      if (!solved && !selected) {
        selected = d;
        wordCircleG.classed('selected', false);
        d3.select(this).classed('selected', true);
      }
    })
    .on('drag', function(event, d) {
      if (!d._didDrag) {
        const dx = (event.x - (d._dragStartX ?? event.x));
        const dy = (event.y - (d._dragStartY ?? event.y));
        if (Math.hypot(dx, dy) < DRAG_THRESHOLD) return; // still a click, don't start dragging yet
        d._didDrag = true;
        if (!event.active) simulation.alphaTarget(BACKGROUND_ALPHA).restart();
        d.fx = d.x; d.fy = d.y;
      }
      d.fx = event.x; d.fy = event.y;
    })
    .on('end', function(event, d) {
      if (d._didDrag) {
        d.fx = null; d.fy = null;
        if (!event.active) simulation.alphaTarget(BACKGROUND_ALPHA);
        // Prevent the subsequent click only if an actual drag occurred
        try { event.sourceEvent && event.sourceEvent.preventDefault && event.sourceEvent.preventDefault(); } catch {}
        snapshotAssignmentOnly();
      }
      d._dragStartX = d._dragStartY = undefined;
      d._didDrag = false;
    });
  wordCircleG.call(dragWordOverlay);
  // Ensure text inside circles doesn't intercept pointer events
  try { wordCircleG.selectAll('text').attr('pointer-events', 'none'); } catch {}

  // Maintain mapping from node alias to assigned word
  const nodeAliasToWord = new Map(assignment);
  const wordCircleByAlias = new Map();
  wordCircleG.each(function(d) { wordCircleByAlias.set(d.alias, d3.select(this)); });
  // Persist initial assignment
  if (persist && typeof persist.save === 'function') {
    const asn0 = {}; nodeAliasToWord.forEach((v,k)=>asn0[k]=v);
    persist.save({ assignment: asn0, solved });
  }

  // Build helper maps for highlighting
  const categoryNodeByAlias = new Map();
  catNodes.each(function(d){ categoryNodeByAlias.set(d.alias, d3.select(this)); });
  const linksByCategory = new Map();
  link.each(function(l){
    const cat = typeof l.target === 'string' ? l.target : l.target.alias;
    const arr = linksByCategory.get(cat) || [];
    arr.push(d3.select(this));
    linksByCategory.set(cat, arr);
  });

  function commonCategoryForAlias(catAlias) {
    // Compute the intersection of categories among assigned words connected to this category node
    const connectedWordAliases = links
      .filter(l => (typeof l.target === 'string' ? l.target : l.target.alias) === catAlias)
      .map(l => (typeof l.source === 'string' ? l.source : l.source.alias));
    if (connectedWordAliases.length === 0) return null;
    let intersection = null;
    for (const wa of connectedWordAliases) {
      const word = nodeAliasToWord.get(wa);
      const catsArr = wordToCategories[word] || [];
      const cats = new Set(catsArr);
      if (!intersection) {
        intersection = cats;
      } else {
        const next = new Set();
        for (const c of intersection) if (cats.has(c)) next.add(c);
        intersection = next;
      }
      if (!intersection || intersection.size === 0) return null;
    }
    // Restrict to categories present in this puzzle only
    const filtered = new Set();
    for (const c of intersection) if (puzzleCategories.has(c)) filtered.add(c);
    if (filtered.size === 0) return null;
    // Prefer a category that has a known emoji, otherwise any
    let chosen = null;
    for (const c of filtered) { if (categoryEmojis[c]) { chosen = c; break; } }
    if (!chosen) chosen = filtered.values().next().value || null;
    return chosen;
  }

  // Helper: build and open the Share modal with current data; also log solve order
  function openShareModalWithData() {
    try {
      console.log('[Nodeword] Opening share modal (helper)');
      // Show modal immediately; then populate
      const shareModal = document.getElementById('shareModal');
      if (!shareModal) { console.warn('[Nodeword] Share modal element not found'); return; }
      shareModal.hidden = false; console.log('[Nodeword] Share modal shown');
      const dNow = new Date();
      const today = `${dNow.getFullYear()}-${String(dNow.getMonth()+1).padStart(2,'0')}-${String(dNow.getDate()).padStart(2,'0')}`;
      // Highlighted word is W1's canonical id
      const w1 = aliasGraph.nodes.find(n=>n.type==='word' && n.alias==='W1');
      const highlightedWord = w1 ? w1.id : '';
      // Persist this level emoji order if missing and build chain
      try {
        const raw = localStorage.getItem('nodeword_state_v1');
        let s = {};
        try { s = raw ? JSON.parse(raw) : {}; console.log('[Nodeword] Share state loaded'); } catch (e) { console.warn('[Nodeword] Bad saved state JSON; ignoring'); s = {}; }
        if (!s.dailyEmojis) s.dailyEmojis = {};
        const level = Math.min(5, Math.max(1, Number(persist && persist.meta && persist.meta.levelToday || 1)));
        const levelKey = 'L' + level;
        if (Array.isArray(solveOrderEmojis) && solveOrderEmojis.length && (!Array.isArray(s.dailyEmojis[levelKey]) || s.dailyEmojis[levelKey].length === 0)) {
          s.dailyEmojis[levelKey] = solveOrderEmojis.slice();
          try { s.version = STORAGE_VERSION; } catch {}
          // Preserve existing dailyEmojis via writeState merge to avoid wiping nuclear-saved chains
          writeState(s);
          // Mirror dailyEmojis into in-memory state so later writes (e.g., Next) don't wipe them
          try { appState.dailyEmojis = s.dailyEmojis; } catch {}
        }
        // Build chain from dailyEmojis; fallback to saved puzzle's solveOrder (by cats → emojis) if empty on L5
        const parts = [];
        for (let i = 1; i <= 5; i++) {
          const arr = s.dailyEmojis['L' + i];
          if (Array.isArray(arr) && arr.length) parts.push(arr.join(''));
        }
        let chain = parts.join('→');
        if (!chain && level === 5) {
          try {
            const saved = (typeof appState === 'object' && appState.puzzle) ? appState.puzzle : null;
            if (saved) {
              if (Array.isArray(saved.solveOrderCats) && saved.solveOrderCats.length) {
                // Map category aliases to emojis
                const aliasToEmoji = new Map();
                for (const node of aliasGraph.nodes) {
                  if (node.type === 'category') aliasToEmoji.set(node.id, categoryEmojis[node.id] || '');
                }
                const ems = saved.solveOrderCats.map(ca => aliasToEmoji.get(ca) || '').filter(Boolean);
                if (ems.length) chain = ems.join('');
              }
              if (!chain && Array.isArray(saved.solveOrderEmojis) && saved.solveOrderEmojis.length) {
                chain = saved.solveOrderEmojis.join('');
              }
            }
          } catch {}
        }
        if (!chain && Array.isArray(solveOrderEmojis) && solveOrderEmojis.length) {
          chain = solveOrderEmojis.join('');
        }
        console.log('[Nodeword] Share chain', { chain: chain || '(empty)', L1: s.dailyEmojis.L1||[], L2: s.dailyEmojis.L2||[], L3: s.dailyEmojis.L3||[], L4: s.dailyEmojis.L4||[], L5: s.dailyEmojis.L5||[] });
        const includeLinkEl = document.getElementById('includeLink');
        const shareTextEl = document.getElementById('shareText');
        console.log('[Nodeword] Share DOM presence', { shareText: !!shareTextEl, includeLink: !!includeLinkEl });
        function buildShare(includeLink) {
          const base = `Nodeword ${today}: "${highlightedWord}"\n\n${chain}\n`;
          return includeLink ? `${base}\nhttps://nodeword.com` : base;
        }
        try {
          if (shareTextEl) { shareTextEl.textContent = buildShare(includeLinkEl ? includeLinkEl.checked : true); console.log('[Nodeword] Share text populated', shareTextEl.textContent); }
          if (includeLinkEl) includeLinkEl.onchange = () => { if (shareTextEl) shareTextEl.textContent = buildShare(includeLinkEl.checked); };
          if (!shareTextEl) console.warn('[Nodeword] shareText element missing');
        } catch (e) { console.warn('[Nodeword] Error wiring share text', e); }
        const copyBtn = document.getElementById('copyShare');
        if (copyBtn) {
          copyBtn.onclick = async () => {
            try { await navigator.clipboard.writeText((shareTextEl && shareTextEl.textContent) || ''); copyBtn.textContent = 'Copied!'; setTimeout(()=>{copyBtn.textContent='Copy';}, 1200);} catch {}
          };
        }
        const closeShare = document.getElementById('closeShare');
        if (closeShare) { closeShare.onclick = ()=>{ shareModal.hidden = true; }; }
      } catch (e) { console.warn('[Nodeword] Share modal population failed', e); }
    } catch {}
  }
  try { window.__nodewordOpenShare = openShareModalWithData; } catch {}

  // Track per-puzzle category solve order (by first time a category becomes valid)
  const solvedCatAliases = new Set();
  // Prefer restoring category aliases; fallback to legacy emojis (will be mapped later if needed)
  let solveOrderCatAliases = (persist && persist.restore && Array.isArray(persist.restore.solveOrderCats)) ? persist.restore.solveOrderCats.slice() : [];
  let solveOrderEmojis = (persist && persist.restore && Array.isArray(persist.restore.solveOrderEmojis)) ? persist.restore.solveOrderEmojis.slice() : [];

  function updateCategoryHighlights() {
    for (const [catAlias, catSel] of categoryNodeByAlias.entries()) {
      const common = commonCategoryForAlias(catAlias);
      const highlight = Boolean(common);
      catSel.classed('highlight', highlight);
      const emojiChar = common ? (categoryEmojis[common] || '✅') : '';
      const textSel = catSel.select('text.cat-emoji');
      textSel.text(emojiChar);
      // Category name label only after entire puzzle is solved
      const labelSel = catSel.select('text.cat-label');
      labelSel.text(highlight && solved ? (common || '') : '');
      const lines = linksByCategory.get(catAlias) || [];
      for (const ln of lines) ln.classed('highlight', highlight);

      // Record first-time solve order for this category (by alias, not emoji)
      if (highlight && !solvedCatAliases.has(catAlias)) {
        solvedCatAliases.add(catAlias);
        solveOrderCatAliases.push(catAlias);
        // Maintain legacy emoji list for backward compatibility if desired
        if (emojiChar) solveOrderEmojis.push(emojiChar); else solveOrderEmojis.push('✅');
        if (persist && typeof persist.save === 'function') {
          try {
            console.log('[Nodeword] solveOrder updated', { cats: solveOrderCatAliases, emojis: solveOrderEmojis });
            persist.save({ solveOrderCats: solveOrderCatAliases, solveOrderEmojis });
          } catch {}
        }
      }
    }
  }

  wordCircleG.on('click', function(event, d) {
    event.stopPropagation();
    // If this click follows a drag, ignore to prevent select-then-deselect flicker
    if (event.defaultPrevented) return;
    const me = d3.select(this);
    if (!selected) {
      selected = d;
      me.classed('selected', true);
      return;
    }
    if (selected && selected.alias === d.alias) {
      // Keep it selected; do not toggle off on second click
      return;
    }
    // Swap words between selected node and current node
    const a = selected.alias;
    const b = d.alias;
    const wa = nodeAliasToWord.get(a);
    const wb = nodeAliasToWord.get(b);
    nodeAliasToWord.set(a, wb);
    nodeAliasToWord.set(b, wa);

    // Re-measure circles for both if needed
    const ma = measureWordCircle(wb);
    const mb = measureWordCircle(wa);
    circleMetrics.set(a, ma);
    circleMetrics.set(b, mb);

    // Update visuals for both circles (text + radius)
    const ga = wordCircleByAlias.get(a);
    const gb = wordCircleByAlias.get(b);
    if (ga) {
      ga.select('circle').transition().duration(160).attr('r', ma.radius);
      ga.selectAll('text').remove();
      for (let i = 0; i < ma.parts.length; i++) {
        ga.append('text').attr('y', ma.startY + i * (ma.fontSize + 2)).style('font-size', ma.fontSize + 'px').text(ma.parts[i]);
      }
    }
    if (gb) {
      gb.select('circle').transition().duration(160).attr('r', mb.radius);
      gb.selectAll('text').remove();
      for (let i = 0; i < mb.parts.length; i++) {
        gb.append('text').attr('y', mb.startY + i * (mb.fontSize + 2)).style('font-size', mb.fontSize + 'px').text(mb.parts[i]);
      }
    }

    // Nudge simulation to settle with new radii
    simulation.alpha(0.18).restart();

    // Visual: selected state cleared, and ensure both circles fly to their node positions on next tick
    deselect();
    // Update highlighted word to follow the first puzzle word on final level
    try {
      const isFinalLevel = Boolean(persist && persist.meta && persist.meta.levelToday === 5);
      if (isFinalLevel) {
        const firstWord = aliasGraph.nodes.filter(n=>n.type==='word')[0]?.id;
        if (firstWord) {
          let targetAlias = null;
          nodeAliasToWord.forEach((w, alias) => { if (w === firstWord) targetAlias = alias; });
          if (targetAlias) {
            wordCircleG.classed('highlighted', false);
            wordCircleG.filter(d => d.alias === targetAlias).classed('highlighted', true);
          }
        }
      }
    } catch {}
    updateCategoryHighlights();
    snapshotAssignmentOnly();
    checkForSolved();
  });

  function setInteractivity(enabled) {
    wordCircleG.style('pointer-events', enabled ? 'all' : 'none');
  }

  // Click anywhere else to clear selection
  svg.on('click', () => { if (!solved) deselect(); });
  // Initial check
  updateCategoryHighlights();
    if (solved) {
    setInteractivity(false);
    const statusEl = document.getElementById('status');
    if (statusEl) statusEl.textContent = 'Puzzle solved!';
    const btn = document.getElementById('nextBtn');
    const shareBtn = document.getElementById('shareBtn');
    // Determine final level (5) using meta passed from caller
    const levelFromTarget = Math.max(1, Number(persist && persist.meta && persist.meta.levelToday || 1));
    const isFinalLevel = levelFromTarget >= 5;
    console.log('[Nodeword] Solved (initial). level=', levelFromTarget, 'final=', isFinalLevel);
      if (isFinalLevel) {
      console.log('[Nodeword] Trigger share modal (initial solved)');
      if (btn) btn.style.display = 'none';
      if (shareBtn) shareBtn.style.display = 'inline-block';
      openShareModalWithData();
    } else {
      console.log('[Nodeword] Non-final level on initial solved; showing Next');
      if (btn) btn.style.display = 'inline-block';
    }
  } else {
    checkForSolved();
  }

  const computePadding = () => Math.max(8, Math.round(Math.min(width, height) * 0.03));
  let basePadding = computePadding();
  const boundaryK = 0.06; // gentler boundary pushback

  function boundaryForce() {
    // Custom force to nudge nodes back within bounds
    for (const node of nodes) {
      const r = radius(node);
      // Inflate boundary considering overlay circle size for word nodes
      const extra = node.type === 'word' ? (circleMetrics.get(node.alias)?.radius || 0) : 0;
      if (node.x < basePadding + r + extra) node.vx += (basePadding + r + extra - node.x) * boundaryK;
      if (node.x > width - basePadding - r - extra) node.vx += (width - basePadding - r - extra - node.x) * boundaryK;
      if (node.y < basePadding + r + extra) node.vy += (basePadding + r + extra - node.y) * boundaryK;
      if (node.y > height - basePadding - r - extra) node.vy += (height - basePadding - r - extra - node.y) * boundaryK;
    }
  }

  const simulation = d3.forceSimulation(nodes)
    .velocityDecay(0.85)
    .alphaDecay(0.08)
    .alphaMin(0.02)
    .force('link', d3.forceLink(links).id(d => d.alias).distance(104).strength(0.6))
    .force('charge', d3.forceManyBody().strength(d => d.type === 'word' ? -220 - (circleMetrics.get(d.alias)?.radius || 0) * 1.2 : -180))
    .force('collide', d3.forceCollide().radius(d => {
      const extra = d.type === 'word' ? (circleMetrics.get(d.alias)?.radius || 0) : 0;
      return radius(d) + 12 + extra;
    }).strength(0.85))
    // Light lane separation: words to left band, categories to right band
    .force('laneX', d3.forceX(d => d.type === 'word' ? (width * 0.35) : (width * 0.65)).strength(0.03))
    .force('catCentroid', categoryCentroidForce(0.03))
    .force('antiCross', antiCrossingForce(0.06))
    .force('boundary', boundaryForce);

  // Keep a gentle background tick so forces continue to act
  const BACKGROUND_ALPHA = 0.02;
  simulation.alphaTarget(BACKGROUND_ALPHA);

  // Nuclear option: iterative auto-fit scaling
  let lastScale = 1;
  const scaleFloor = 0.55;
  const scaleCeil = 1.05;
  const scaleStep = 0.008; // gentler per-tick scaling
  const scaleEps = 0.02; // wider hysteresis to avoid jitter
  function fitSceneIfNeeded() {
    let minX = Infinity, minY = Infinity, maxX = -Infinity, maxY = -Infinity;
    for (const n of nodes) {
      const r = radius(n);
      const extra = n.type === 'word' ? (circleMetrics.get(n.alias)?.radius || 0) : 0;
      const x1 = (n.x || 0) - (r + extra);
      const y1 = (n.y || 0) - (r + extra);
      const x2 = (n.x || 0) + (r + extra);
      const y2 = (n.y || 0) + (r + extra);
      if (x1 < minX) minX = x1;
      if (y1 < minY) minY = y1;
      if (x2 > maxX) maxX = x2;
      if (y2 > maxY) maxY = y2;
    }
    const contentW = Math.max(1, maxX - minX);
    const contentH = Math.max(1, maxY - minY);
    const sx = (width - 2 * basePadding) / contentW;
    const sy = (height - 2 * basePadding) / contentH;
    const target = Math.min(sx, sy);
    let next = lastScale;
    if (target < lastScale - scaleEps) {
      next = Math.max(scaleFloor, lastScale - scaleStep);
    } else if (target > lastScale + scaleEps) {
      next = Math.min(scaleCeil, lastScale + scaleStep);
    }
    if (Math.abs(next - lastScale) > 1e-3 || lastScale !== 1) {
      const cx = (minX + maxX) / 2;
      const cy = (minY + maxY) / 2;
      // Bias scene slightly left by 4px to compensate for subpixel rounding on narrow screens
      scene.attr('transform', `translate(${width / 2 - 4},${height / 2}) scale(${next}) translate(${-cx},${-cy})`);
      lastScale = next;
    }
  }

  let firstTickLogged = false;
  simulation.alpha(0); // start paused to let KK seed positions be visible
  setTimeout(() => simulation.alpha(0.4).restart(), 120); // then relax
  // Remove periodic snapshots; we'll persist on explicit events only
  function snapshotAssignmentOnly() {
    try {
      if (!persist || typeof persist.save !== 'function') return;
      const asn = {}; nodeAliasToWord.forEach((v,k)=>asn[k]=v);
      persist.save({ assignment: asn });
    } catch {}
  }

  simulation.on('tick', () => {
    link
      .attr('x1', d => d.source.x)
      .attr('y1', d => d.source.y)
      .attr('x2', d => d.target.x)
      .attr('y2', d => d.target.y);
    // Clamp final positions to stay inside the viewport
    node
      .attr('transform', d => {
        const r = radius(d);
        const extra = d.type === 'word' ? (circleMetrics.get(d.alias)?.radius || 0) : 0;
        d.x = Math.max(basePadding + r + extra, Math.min(width - basePadding - r - extra, d.x));
        d.y = Math.max(basePadding + r + extra, Math.min(height - basePadding - r - extra, d.y));
        return `translate(${d.x},${d.y})`;
      });

    // Keep word circles centered atop their assigned word nodes
    wordCircleG
      .attr('transform', d => `translate(${d.x},${d.y})`);
    fitSceneIfNeeded();
    if (!firstTickLogged) {
      firstTickLogged = true;
      console.log('[Nodeword] Force layout started ticking');
    }
  });

  // Resize handling for responsiveness
  const onResize = () => {
    const prevW = width;
    const prevH = height;
    width = container.clientWidth;
    height = Math.min(820, Math.max(520, Math.floor(window.innerHeight * 0.75)));
    basePadding = computePadding();
    svg.attr('width', width).attr('height', height);
    // Update lane force target on resize
    simulation.force('laneX', d3.forceX(d => d.type === 'word' ? (width * 0.35) : (width * 0.65)).strength(0.03));
    // Immediately translate nodes toward the new center so they don't sit offscreen
    const dx = (width - prevW) / 2;
    const dy = (height - prevH) / 2;
    if (dx !== 0 || dy !== 0) {
      for (const n of nodes) {
        n.x = (n.x || 0) + dx;
        n.y = (n.y || 0) + dy;
        if (n.fx != null) n.fx += dx;
        if (n.fy != null) n.fy += dy;
      }
    }
    // Reset scene scale and refit after resize
    scene.attr('transform', null);
    lastScale = 1;
    simulation.alpha(0.8).restart();
  };
  window.addEventListener('resize', onResize, { passive: true });

  function checkForSolved() {
    // Count diamonds without a valid common category under current assignment
    let unsolvedCount = 0;
    for (const catAlias of categoryNodeByAlias.keys()) {
      if (!commonCategoryForAlias(catAlias)) unsolvedCount++;
    }
    if (unsolvedCount === 0) {
      if (!solved) {
        console.log('[Nodeword] SOLVED: all diamonds complete');
        solved = true;
        setInteractivity(false);
        const statusEl = document.getElementById('status');
        if (statusEl) statusEl.textContent = 'Puzzle solved!';
        const btn = document.getElementById('nextBtn');
        const shareBtn = document.getElementById('shareBtn');
        // Mark solved via persistence adapter for unified state flow
        try {
          if (persist && typeof persist.save === 'function')
            persist.save({ solved: true, solveOrderCats: (solveOrderCatAliases||[]).slice() });
        } catch {}
        // If daily cap reached, show Share instead of Next and open share modal
        const dailyLevelCap = 5;
        let reachedCap = false;
        // Reveal labels and ensure emoji/text reflects the chosen intersection category
        updateCategoryHighlights();
        // Increment stats immediately upon solve
        try {
          const raw = localStorage.getItem('nodeword_state_v1');
          const s = raw ? JSON.parse(raw) : {};
          const today = (()=>{const d=new Date();return `${d.getFullYear()}-${d.getMonth()+1}-${d.getDate()}`;})();
          if (s.lastDay !== today) { s.today = 0; s.lastDay = today; }
          s.total = (s.total||0)+1; s.today=(s.today||0)+1;
          // perfect day tracking
          if (s.today === dailyLevelCap) {
            s.perfectDays = (s.perfectDays||0) + 1;
            s.streak = (s.streak||0) + 1;
            s.completedPerfectDay = true;
            reachedCap = true;
          }
          // Ensure puzzle snapshot exists even if earlier saves were gated by consent
          try {
            const allowPersist = (typeof appState === 'object') ? (appState.consent !== false) : true;
            if (allowPersist) {
              const wordsList = aliasGraph.nodes.filter(n=>n.type==='word').map(n=>n.id);
              const catList = aliasGraph.nodes.filter(n=>n.type==='category').map(n=>n.id);
              // Map alias->id for words and categories
              const aliasToId = new Map(aliasGraph.nodes.map(n=>[n.alias, n.id]));
              const edgesList = aliasGraph.links.map(l=>({
                word: aliasToId.get(typeof l.source==='string'? l.source : l.source.alias),
                category: aliasToId.get(typeof l.target==='string'? l.target : l.target.alias)
              })).filter(e=>e.word && e.category);
              const snapshot = { target: Math.min(Math.max(5, currentTargetWords), maxTargetWords), graph: { words: wordsList, categories: catList, edges: edgesList } };
              const asnObj = {}; nodeAliasToWord.forEach((v,k)=>asnObj[k]=v);
              // Persist per-level solve order by category aliases and mapped emojis into puzzle snapshot
              const aliasToEmoji = new Map();
              for (const node of aliasGraph.nodes) if (node.type==='category') aliasToEmoji.set(node.alias, categoryEmojis[node.id]||'');
              const thisLevelCats = (solveOrderCatAliases||[]).slice();
              const thisLevelEmojis = thisLevelCats.map(ca=>aliasToEmoji.get(ca)||'').filter(Boolean);
              s.puzzle = { ...(s.puzzle||{}), ...snapshot, assignment: asnObj, solved:true, solveOrderCats: thisLevelCats, solveOrderEmojis: thisLevelEmojis };
              try {
                console.log('[Nodeword] Saving solved state to storage', {
                  target: snapshot.target,
                  words: snapshot.graph.words.length,
                  categories: snapshot.graph.categories.length,
                  edges: snapshot.graph.edges.length
                });
              } catch {}
            }
          } catch {}
          // Nuclear option: Always persist dailyEmojis for this level outside consent gate using a deep-clone write
          try {
            const aliasToEmoji2 = new Map();
            for (const node of aliasGraph.nodes) if (node.type==='category') aliasToEmoji2.set(node.alias, categoryEmojis[node.id]||'');
            const catsNow = (solveOrderCatAliases||[]).slice();
            const emsNow = catsNow.map(ca=>aliasToEmoji2.get(ca)||'').filter(Boolean);
            // Derive level from current graph word count to avoid relying on outer variables
            const wordsCountNow = aliasGraph.nodes.filter(n=>n.type==='word').length;
            const levelTodayNow = Math.min(5, Math.max(1, (wordsCountNow - 4)));
            const prevRaw = localStorage.getItem('nodeword_state_v1');
            const prev = prevRaw ? JSON.parse(prevRaw) : {};
            if (!prev.dailyEmojis) prev.dailyEmojis = {};
            prev.dailyEmojis['L'+levelTodayNow] = (emsNow || []).slice();
            // Deep-clone snapshot to avoid pointer aliasing
            localStorage.setItem('nodeword_state_v1', JSON.stringify(JSON.parse(JSON.stringify(prev))));
            console.log('[Nodeword] Nuclear save: dailyEmojis updated', { level: levelTodayNow, emojis: prev.dailyEmojis['L'+levelTodayNow], cats: catsNow });
          } catch (e) { console.warn('[Nodeword] Nuclear save failed', e); }
          try { s.version = STORAGE_VERSION; } catch {}
          writeState(s);
          try { logDailyDebug('post-solve'); } catch {}
          // Keep in-memory state aligned so later persistence (e.g., node position snapshots) does not overwrite solved flag
          try {
            appState.puzzle = { ...(appState.puzzle||{}), ...snapshot, assignment: asnObj, solved: true };
            appState.target = snapshot.target;
            // Persist updated appState (including mirrored dailyEmojis)
            writeState(appState);
            console.log('[Nodeword] In-memory state updated to solved for current target', {
              target: appState.target,
              hasAssignment: !!appState.puzzle.assignment,
              solved: !!appState.puzzle.solved
            });
          } catch {}
          const statTotal = document.getElementById('statTotal'); if (statTotal) statTotal.textContent = String(s.total||0);
          const statStreak = document.getElementById('statStreak'); if (statStreak) statStreak.textContent = String(s.streak||0);
          const statPerfectDays = document.getElementById('statPerfectDays'); if (statPerfectDays) statPerfectDays.textContent = String(s.perfectDays||0);
          // Keep in-memory appState in sync; level stays based on current target until Next is clicked
          try { appState.total = s.total; appState.today = s.today; appState.streak = s.streak; appState.perfectDays = s.perfectDays; } catch {}
          const lvlEl = document.getElementById('levelIndicator');
          if (lvlEl) {
            const levelFromTarget = Math.min(5, Math.max(1, ((appState.target || 5) - 4)));
            lvlEl.textContent = `Level ${levelFromTarget}`;
          }
        } catch {}
        // Determine level from target words (final level is 5)
        const levelFromTarget = Math.min(5, Math.max(1, Number(persist && persist.meta && persist.meta.levelToday || 1)));
        const isFinalLevel = levelFromTarget === 5;
        console.log('[Nodeword] SOLVED handler. level=', levelFromTarget, 'final=', isFinalLevel);
        // Toggle buttons based on final level
        if (btn) btn.style.display = isFinalLevel ? 'none' : 'inline-block';
        if (shareBtn) shareBtn.style.display = isFinalLevel ? 'inline-block' : 'none';
        // If final level, open share modal with centralized helper
        if (isFinalLevel) {
          console.log('[Nodeword] Triggering share modal (solve)');
          openShareModalWithData();
        }
      }
    } else {
      console.log('[Nodeword] NOT SOLVED: remaining diamonds', unsolvedCount);
    }
  }
}

document.addEventListener("DOMContentLoaded", () => {
  console.log('[Nodeword] DOMContentLoaded');
  if (typeof d3 === 'undefined') {
    console.error('[Nodeword] D3 not loaded. Force graph will not render.');
  } else {
    console.log('[Nodeword] D3 version', d3.version);
  }
  const title = document.querySelector("h1");
  if (title) {
    title.animate(
      [ { transform: "translateY(8px)", opacity: 0 }, { transform: "translateY(0)", opacity: 1 } ],
      { duration: 600, easing: "cubic-bezier(.2,.8,.2,1)", fill: "both" }
    );
  }

  const container = document.getElementById('puzzle');
  const nextBtn = document.getElementById('nextBtn');
  const statsBtn = document.getElementById('statsBtn');
  const shareBtn = document.getElementById('shareBtn');
  const statsModal = document.getElementById('statsModal');
  const closeStats = document.getElementById('closeStats');
  const statTotal = document.getElementById('statTotal');
  const statToday = document.getElementById('statToday');
  const statBest = document.getElementById('statBest');
  const consent = document.getElementById('consent');
  const acceptConsent = document.getElementById('acceptConsent');
  const declineConsent = document.getElementById('declineConsent');
  // Persistence helpers and stats
  const STORAGE_KEY = 'nodeword_state_v1';
  const STORAGE_VERSION = 2;
  function readState() {
    try { const raw = localStorage.getItem(STORAGE_KEY); return raw ? JSON.parse(raw) : null; } catch { return null; }
  }
  function canPersist() { return appState && appState.consent === true; }
  function writeState(s) {
    if (!canPersist()) return;
    try {
      // Preserve/merge previously saved dailyEmojis if caller didn't include them or provided empties
      try {
        const prevRaw = localStorage.getItem(STORAGE_KEY);
        const prev = prevRaw ? JSON.parse(prevRaw) : null;
        if (prev && prev.dailyEmojis) {
          if (!s.dailyEmojis) s.dailyEmojis = {};
          const keys = ['L1','L2','L3','L4','L5'];
          for (const k of keys) {
            const cur = s.dailyEmojis[k];
            const old = prev.dailyEmojis[k];
            if (!Array.isArray(cur) || cur.length === 0) {
              if (Array.isArray(old) && old.length > 0) s.dailyEmojis[k] = old.slice();
            }
          }
        }
      } catch {}
      s.version = STORAGE_VERSION;
      localStorage.setItem(STORAGE_KEY, JSON.stringify(s));
      try {
        const dbg = (s && s.dailyEmojis) ? { L1: s.dailyEmojis.L1||[], L2: s.dailyEmojis.L2||[], L3: s.dailyEmojis.L3||[], L4: s.dailyEmojis.L4||[], L5: s.dailyEmojis.L5||[] } : null;
        console.log('[Nodeword] writeState saved dailyEmojis', dbg);
      } catch {}
    } catch {}
  }
  function pad2(n){ return String(n).padStart(2,'0'); }
  function todayStamp() {
    const d = new Date();
    return `${d.getFullYear()}-${pad2(d.getMonth()+1)}-${pad2(d.getDate())}`;
  }
  function initState() {
    const today = todayStamp();
    const s = readState();
    if (!s) return { version: STORAGE_VERSION, consent: null, total: 0, today: 0, best: 0, lastDay: today, target: 5, puzzle: null, streak: 0, perfectDays: 0, dailyEmojis: {} };
    // Migration: if missing/older version or clearly invalid shape, wipe state
    try {
      const invalid = (
        !('version' in s) || (typeof s.version !== 'number') || s.version < STORAGE_VERSION ||
        (s.puzzle && (!s.puzzle.graph || !Array.isArray(s.puzzle.graph?.words) || !Array.isArray(s.puzzle.graph?.categories)))
      );
      if (invalid) {
        try { localStorage.removeItem(STORAGE_KEY); } catch {}
        return { version: STORAGE_VERSION, consent: null, total: 0, today: 0, best: 0, lastDay: today, target: 5, puzzle: null, streak: 0, perfectDays: 0, dailyEmojis: {} };
      }
    } catch {}
    if (s.lastDay !== today) {
      // Reset daily counters; keep streak only if yesterday completed all five
      if (!s.completedPerfectDay) s.streak = 0; // break streak if not perfect yesterday
      s.completedPerfectDay = false;
      s.today = 0; s.lastDay = today; s.target = 5; s.puzzle = null; s.dailyEmojis = {};
    }
    // Ensure version field is stamped going forward
    if (s.version !== STORAGE_VERSION) s.version = STORAGE_VERSION;
    return s;
  }
  let appState = initState();
  function logDailyDebug(tag) {
    try {
      const raw = localStorage.getItem(STORAGE_KEY);
      const s = raw ? JSON.parse(raw) : null;
      if (s && s.dailyEmojis) {
        const arrays = { L1: s.dailyEmojis.L1||[], L2: s.dailyEmojis.L2||[], L3: s.dailyEmojis.L3||[], L4: s.dailyEmojis.L4||[], L5: s.dailyEmojis.L5||[] };
        const chains = {
          L1: Array.isArray(arrays.L1) ? arrays.L1.join('') : '',
          L2: Array.isArray(arrays.L2) ? arrays.L2.join('') : '',
          L3: Array.isArray(arrays.L3) ? arrays.L3.join('') : '',
          L4: Array.isArray(arrays.L4) ? arrays.L4.join('') : '',
          L5: Array.isArray(arrays.L5) ? arrays.L5.join('') : ''
        };
        console.log(`[Nodeword] Daily emoji arrays ${tag? '('+tag+')':''}`, arrays);
        console.log(`[Nodeword] Daily emoji chains ${tag? '('+tag+')':''}`, chains);
      } else {
        console.log('[Nodeword] No dailyEmojis found in storage', { tag });
      }
    } catch (e) { console.warn('[Nodeword] logDailyDebug failed', e); }
  }
  try {
    console.log('[Nodeword] appState initialized', {
      target: appState.target,
      today: appState.today,
      total: appState.total,
      consent: appState.consent,
      puzzle: appState.puzzle ? { solved: !!appState.puzzle.solved, target: appState.puzzle.target } : null
    });
    // Log any saved daily emoji chains and saved solveOrderEmojis for debugging
    try {
      const rawDbg = localStorage.getItem(STORAGE_KEY);
      let sDbg = null;
      try { sDbg = rawDbg ? JSON.parse(rawDbg) : null; } catch {}
      if (sDbg && sDbg.dailyEmojis) {
        console.log('[Nodeword] Daily emoji arrays', {
          L1: sDbg.dailyEmojis.L1 || [],
          L2: sDbg.dailyEmojis.L2 || [],
          L3: sDbg.dailyEmojis.L3 || [],
          L4: sDbg.dailyEmojis.L4 || [],
          L5: sDbg.dailyEmojis.L5 || []
        });
        try {
          const chains = {
            L1: Array.isArray(sDbg.dailyEmojis.L1) ? sDbg.dailyEmojis.L1.join('') : '',
            L2: Array.isArray(sDbg.dailyEmojis.L2) ? sDbg.dailyEmojis.L2.join('') : '',
            L3: Array.isArray(sDbg.dailyEmojis.L3) ? sDbg.dailyEmojis.L3.join('') : '',
            L4: Array.isArray(sDbg.dailyEmojis.L4) ? sDbg.dailyEmojis.L4.join('') : '',
            L5: Array.isArray(sDbg.dailyEmojis.L5) ? sDbg.dailyEmojis.L5.join('') : ''
          };
          console.log('[Nodeword] Daily emoji chains', chains);
        } catch {}
      } else {
        console.log('[Nodeword] No dailyEmojis found in storage');
      }
      if (appState.puzzle && Array.isArray(appState.puzzle.solveOrderEmojis)) {
        console.log('[Nodeword] Saved puzzle solveOrderEmojis', appState.puzzle.solveOrderEmojis);
      } else {
        console.log('[Nodeword] No saved puzzle solveOrderEmojis');
      }
    } catch {}
    logDailyDebug('on-load');
  } catch {}
  function updateStatsUI() {
    if (statTotal) statTotal.textContent = String(appState.total || 0);
    const statStreak = document.getElementById('statStreak'); if (statStreak) statStreak.textContent = String(appState.streak || 0);
    const statPerfectDays = document.getElementById('statPerfectDays'); if (statPerfectDays) statPerfectDays.textContent = String(appState.perfectDays || 0);
    const lvlEl = document.getElementById('levelIndicator');
    const levelFromTarget = Math.min(5, Math.max(1, ((appState.target || 5) - 4)));
    if (lvlEl) lvlEl.textContent = `Level ${levelFromTarget}`;
  }
  updateStatsUI();
  function showConsentIfNeeded() { if (consent) consent.hidden = !(appState.consent === null); }
  showConsentIfNeeded();
  acceptConsent?.addEventListener('click', () => { appState.consent = true; writeState(appState); showConsentIfNeeded(); });
  declineConsent?.addEventListener('click', () => { appState.consent = false; try { localStorage.removeItem(STORAGE_KEY); } catch {} showConsentIfNeeded(); });
  function openStats() { if (statsModal) statsModal.hidden = false; }
  function closeStatsModal() { if (statsModal) statsModal.hidden = true; }
  statsBtn?.addEventListener('click', openStats);
  closeStats?.addEventListener('click', closeStatsModal);
  shareBtn?.addEventListener('click', () => { const m = document.getElementById('shareModal'); if (m) m.hidden = false; });
  // Close stats when clicking backdrop
  statsModal?.addEventListener('click', (e) => { if (e.target === statsModal) closeStatsModal(); });
  // Ensure stats modal starts hidden
  if (statsModal) statsModal.hidden = true;

  // Daily cap: allow up to 5 levels per day; deterministic per day+level
  // Progressive puzzle size still scales up to 12, but only first 5 are playable per day
  let currentTargetWords = appState.target || 5;
  const maxTargetWords = 12;
  const dailyLevelCap = 5;

  let wordData = null;

  async function generateAndRender() {
    container.textContent = 'Generating…';
    console.log('[Nodeword] Generating puzzle…');
    // Hide next button until solved and clear status message
    if (nextBtn) nextBtn.style.display = 'none';
    const statusEl = document.getElementById('status');
    if (statusEl) statusEl.textContent = '';
    // Ensure level indicator reflects current target-based level immediately
    try {
      const lvlToday = Math.min(5, Math.max(1, ((appState.target || 5) - 4)));
      const lvlEl = document.getElementById('levelIndicator');
      if (lvlEl) lvlEl.textContent = `Level ${lvlToday}`;
    } catch {}
    try {
      if (!wordData) wordData = await fetchWordData();
      const cfg = NODEWORD_CONFIG;
      const target = Math.min(Math.max(5, currentTargetWords), maxTargetWords);
      // If a solved puzzle is saved for this target, display it as-is (no generation)
      const saved = appState.puzzle;
      console.log('[Nodeword] Checking saved puzzle for solved restore', {
        exists: Boolean(saved),
        solved: Boolean(saved && saved.solved),
        savedTarget: saved && saved.target,
        target,
        hasGraph: Boolean(saved && saved.graph)
      });
      try { if (saved) console.log('[Nodeword] Saved snapshot (debug)', { 
        keys: Object.keys(saved||{}), 
        hasAssignment: !!saved.assignment,
        hasNodePositions: !!saved.nodePositions,
        graphWords: saved.graph && saved.graph.words && saved.graph.words.length,
        graphCats: saved.graph && saved.graph.categories && saved.graph.categories.length,
        graphEdges: saved.graph && saved.graph.edges && saved.graph.edges.length
      }); } catch {}
      const levelTodayCheck = Math.min(5, Math.max(1, ((target || 5) - 4)));
      const solvedByStats = (appState.today || 0) >= levelTodayCheck;
      if (saved && saved.target === target && saved.graph && (saved.solved === true || solvedByStats)) {
        const graph = saved.graph;
        console.log('[Nodeword] Restoring solved puzzle from save', { solvedByFlag: Boolean(saved && saved.solved), solvedByStats });
        const aliasGraph = toAliasGraph(graph);
        const catEmojis = await fetchCategoryEmojis();
        // Derive level from saved target words to avoid relying on appState
        const levelTodaySaved = Math.min(5, Math.max(1, ((saved.target || 5) - 4)));
        // Ensure in-memory target reflects the current puzzle for Next button progression
        try { currentTargetWords = saved.target; appState.target = saved.target; writeState(appState); } catch {}
        const persist = {
          restore: { nodePositions: saved.nodePositions || null, assignment: null, solved: true },
          save(payload) { /* solved: do nothing */ },
          meta: { levelToday: levelTodaySaved }
        };
        console.log('[Nodeword] Restore meta.levelToday =', levelTodaySaved, 'from saved.target =', saved.target);
        renderForceGraph(container, aliasGraph, wordData, catEmojis, persist);
        console.log('[Nodeword] Graph rendered to SVG (solved restore)');
        return;
      }

      // Validate saved puzzle
      function validateSavedPuzzle(saved) {
        if (!saved || !saved.graph) return false;
        const g = saved.graph;
        if (!Array.isArray(g.words) || !Array.isArray(g.categories) || !Array.isArray(g.edges)) return false;
        for (const w of g.words) { if (!wordData[w] || !Array.isArray(wordData[w])) return false; }
        const catsSet = new Set(g.categories);
        for (const e of g.edges) {
          if (!e || typeof e.word !== 'string' || typeof e.category !== 'string') return false;
          if (!wordData[e.word] || !wordData[e.word].includes(e.category)) return false;
          if (!catsSet.has(e.category)) return false;
        }
        if (saved.assignment) {
          const wordsSet = new Set(g.words);
          for (const [, assignedWord] of Object.entries(saved.assignment)) {
            if (!wordsSet.has(assignedWord)) return false;
          }
        }
        return true;
      }
      function isSolvableByCanonical(g) {
        try {
          const alias = toAliasGraph(g);
          const puzzleCategories = new Set(alias.nodes.filter(n=>n.type==='category').map(n=>n.id));
          const catToWordAliases = new Map();
          for (const l of alias.links) {
            const ca = l.target; const wa = l.source;
            const arr = catToWordAliases.get(ca) || []; arr.push(wa); catToWordAliases.set(ca, arr);
          }
          const wordAliasToWord = new Map();
          alias.nodes.filter(n=>n.type==='word').forEach((n, i)=>{ wordAliasToWord.set(n.alias, g.words[i]); });
          for (const [catAlias, wordAliases] of catToWordAliases.entries()) {
            if (!wordAliases || wordAliases.length === 0) return false;
            let intersection = null;
            for (const wa of wordAliases) {
              const w = wordAliasToWord.get(wa);
              const cats = new Set(wordData[w] || []);
              if (!intersection) { intersection = cats; } else {
                const next = new Set();
                for (const c of intersection) if (cats.has(c)) next.add(c);
                intersection = next;
              }
              if (!intersection || intersection.size === 0) return false;
            }
            let ok = false; for (const c of intersection) { if (puzzleCategories.has(c)) { ok = true; break; } }
            if (!ok) return false;
          }
          return true;
        } catch { return false; }
      }

      async function generateGraphWithRetries(targetWords) {
        const multipliers = [1, 2, 3, 4];
        for (const mult of multipliers) {
          try {
            const g = generatePuzzleGraph(wordData, targetWords, cfg.maxDegree, cfg.maxAttempts * mult);
            if (isSolvableByCanonical(g)) return g;
          } catch (e) {
            console.warn('[Nodeword] generation attempt failed (mult=', mult, '):', e.message || e);
          }
        }
        // Last resort: try once without canonical check but catch errors
        try {
          return generatePuzzleGraph(wordData, targetWords, cfg.maxDegree, cfg.maxAttempts * 5);
        } catch (e) {
          console.error('[Nodeword] final generation failed:', e.message || e);
          return null;
        }
      }

      const savedValid = validateSavedPuzzle(appState.puzzle) && appState.puzzle.target === target && isSolvableByCanonical(appState.puzzle.graph);
      let graph;
      if (savedValid) {
        graph = appState.puzzle.graph;
        console.log('[Nodeword] Using savedValid graph for generation path');
      } else {
        // Deterministic seed by day + level (attempt count provides fallbacks)
        const today = todayStamp();
        // Derive level from target words to keep consistent with word count
        const levelTodayGen = Math.min(5, Math.max(1, (target - 4)));
        console.log('[Nodeword] Generating graph with constraints…', { ...cfg, targetWords: target, seed: `${today}-L${levelTodayGen}` });
        graph = await withSeed(`${today}-L${levelTodayGen}`, async () => generateGraphWithRetries(target));
        if (appState.puzzle) { appState.puzzle = null; writeState(appState); }
        if (!graph) {
          // Give up gracefully: inform and try again with a fresh state
          console.warn('[Nodeword] Unable to generate puzzle after retries, retrying fresh…');
          setTimeout(() => generateAndRender(), 50);
          return;
        }
      }
      console.log('[Nodeword] Graph generated:', {
        words: graph.words.length,
        categories: graph.categories.length,
        edges: graph.edges.length
      });
      // Ensure in-memory target reflects the current puzzle for Next button progression
      try { currentTargetWords = target; appState.target = target; writeState(appState); } catch {}
      const aliasGraph = toAliasGraph(graph);
      const catEmojis = await fetchCategoryEmojis();
      // Provide persistence adapter
      const persist = {
        restore: savedValid ? {
          nodePositions: appState.puzzle.nodePositions || null,
          assignment: appState.puzzle.assignment || null,
          solved: appState.puzzle.solved || false,
        } : null,
        save(payload) {
          // Only log when marking solved to avoid console spam from position snapshots
          try { if (payload && payload.solved === true) console.log('[Nodeword] persist.save solved=true'); } catch {}
          // Merge payload with existing puzzle snapshot; preserve fields unless explicitly overridden
          const prev = appState.puzzle || {};
          const nodePositions = (payload && payload.nodePositions) ?? prev.nodePositions ?? null;
          const assignment = (payload && payload.assignment) ?? prev.assignment ?? null;
          // Once solved becomes true, it stays true
          const solvedFlag = (prev.solved === true) || (payload && payload.solved === true);
          if (solvedFlag === true && !(prev.solved === true)) {
            try { console.log('[Nodeword] Persist adapter: solved flag set true; merging payload to state'); } catch {}
          }
          appState.puzzle = {
            ...prev,
            startedAt: prev.startedAt || Date.now(),
            target,
            graph,
            nodePositions,
            assignment,
            ...payload,
            solved: solvedFlag,
          };
          const nextTarget = Math.max(appState.target || 0, target || 0);
          if (nextTarget !== appState.target) {
            try { console.log('[Nodeword] Updating appState.target', { from: appState.target, to: nextTarget }); } catch {}
            appState.target = nextTarget;
          }
          writeState(appState);
        }
      };
      // Pass level meta to renderer for W1 highlight
      const levelTodayMeta = Math.min(5, Math.max(1, (target - 4)));
      persist.meta = { levelToday: levelTodayMeta };
      renderForceGraph(container, aliasGraph, wordData, catEmojis, persist);
      console.log('[Nodeword] Graph rendered to SVG');
    } catch (err) {
      container.textContent = 'Failed to generate puzzle.';
      console.error(err);
    }
  }

  nextBtn?.addEventListener('click', () => {
    console.log('[Nodeword] Next puzzle clicked');
    logDailyDebug('before-next');
    // Enforce daily cap: allow only 5 solves per day
    if ((appState.today || 0) >= dailyLevelCap) {
      const statusEl = document.getElementById('status');
      if (statusEl) statusEl.textContent = 'Come back tomorrow for more Nodeword!';
      nextBtn.style.display = 'none';
      return;
    }
    // Stats increment on solve happens at solve time; advance difficulty and regenerate
    currentTargetWords = Math.min((appState.target || currentTargetWords) + 1, maxTargetWords);
    appState.target = currentTargetWords;
    // Clear current puzzle restore state for next puzzle
    appState.puzzle = null;
    writeState(appState);
    logDailyDebug('after-next');
    generateAndRender();
  });
  if (nextBtn) nextBtn.style.display = 'none';
  generateAndRender();
});

