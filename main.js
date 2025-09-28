function getRandomInt(max) {
  return Math.floor(Math.random() * max);
}

function choice(array) {
  return array[getRandomInt(array.length)];
}

function shuffle(array) {
  for (let i = array.length - 1; i > 0; i--) {
    const j = Math.floor(Math.random() * (i + 1));
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

      const pickExistingWord = Math.random() < 0.5 || wordSet.size === 1;
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
    const restarts = 3;
    for (let rs = 0; rs < restarts; rs++) {
      // Init positions randomly on a circle with a phase offset per restart
      let x = new Array(N), y = new Array(N);
      const phase = Math.random() * 2 * Math.PI;
      const radius0 = Math.min(width, height) * 0.3;
      for (let i = 0; i < N; i++) {
        const a = phase + (2 * Math.PI * i) / N;
        x[i] = width / 2 + Math.cos(a) * radius0;
        y[i] = height / 2 + Math.sin(a) * radius0;
      }
      // Gradient descent on stress
      const iters = 800;
      const step = 0.0015;
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
            const k = K(i, j);
            if (k === 0) continue;
            const f = k * (r - l) / r;
            const fxij = f * dx;
            const fyij = f * dy;
            fx[i] -= fxij; fy[i] -= fyij;
            fx[j] += fxij; fy[j] += fyij;
          }
        }
        let maxDelta = 0;
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
    // Normalize to padded viewport
    let minX = Math.min(...x), maxX = Math.max(...x);
    let minY = Math.min(...y), maxY = Math.max(...y);
    let cw = Math.max(1, maxX - minX), ch = Math.max(1, maxY - minY);
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
    .attr('class', d => `node ${d.type}`)
    .call(d3.drag()
      .on('start', (event, d) => {
        if (!event.active) simulation.alphaTarget(0.08).restart();
        d.fx = d.x; d.fy = d.y;
      })
      .on('drag', (event, d) => {
        d.fx = event.x; d.fy = event.y;
      })
      .on('end', (event, d) => {
        if (!event.active) simulation.alphaTarget(0);
        d.fx = null; d.fy = null;
      })
    );

  const radius = d => d.type === 'category' ? 16 : 12;

  // Draw shapes: word nodes as circles, category nodes as diamonds (rotated squares)
  node.filter(d => d.type === 'word')
    .append('circle')
    .attr('r', radius);

  const catNodes = node.filter(d => d.type === 'category');
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
  if (persist && persist.restore && persist.restore.assignment) {
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

  // Click handling for swapping circles between nodes
  let selected = null; // d (node datum) for selected word-circle
  let solved = Boolean(persist && persist.restore && persist.restore.solved);

  function deselect() {
    selected = null;
    wordCircleG.classed('selected', false);
  }

  function flyTo(gSel, x, y) {
    gSel.transition().duration(280).ease(d3.easeCubicOut)
      .attr('transform', `translate(${x},${y})`);
  }

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
    }
  }

  wordCircleG.on('click', function(event, d) {
    event.stopPropagation();
    const me = d3.select(this);
    if (!selected) {
      selected = d;
      me.classed('selected', true);
      return;
    }
    if (selected && selected.alias === d.alias) {
      deselect();
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
    simulation.alpha(0.4).restart();

    // Visual: selected state cleared, and ensure both circles fly to their node positions on next tick
    deselect();
    updateCategoryHighlights();
    checkForSolved();
  });

  function setInteractivity(enabled) {
    wordCircleG.style('pointer-events', enabled ? 'all' : 'none');
  }

  // Click anywhere else to clear selection
  svg.on('click', () => { if (!solved) deselect(); });
  // Initial check
  updateCategoryHighlights();
  checkForSolved();

  const computePadding = () => Math.max(24, Math.round(Math.min(width, height) * 0.06));
  let basePadding = computePadding();
  const boundaryK = 0.08; // gentler boundary pushback

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
    .force('link', d3.forceLink(links).id(d => d.alias).distance(72).strength(0.75))
    .force('charge', d3.forceManyBody().strength(d => d.type === 'word' ? -140 - (circleMetrics.get(d.alias)?.radius || 0) * 1.2 : -140))
    .force('collide', d3.forceCollide().radius(d => {
      const extra = d.type === 'word' ? (circleMetrics.get(d.alias)?.radius || 0) : 0;
      return radius(d) + 10 + extra;
    }).strength(0.9))
    .force('center', d3.forceCenter(width / 2, height / 2))
    .force('boundary', boundaryForce);

  // Nuclear option: iterative auto-fit scaling
  let lastScale = 1;
  const scaleFloor = 0.55;
  const scaleCeil = 1.05;
  const scaleStep = 0.03; // adjust per tick
  const scaleEps = 0.01; // hysteresis to avoid jitter
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
      scene.attr('transform', `translate(${width / 2},${height / 2}) scale(${next}) translate(${-cx},${-cy})`);
      lastScale = next;
    }
  }

  let firstTickLogged = false;
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
    // Update center force to new viewport
    simulation.force('center', d3.forceCenter(width / 2, height / 2));
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
        if (btn) btn.style.visibility = 'visible';
        // Reveal labels and ensure emoji/text reflects the chosen intersection category
        updateCategoryHighlights();
        // Increment stats immediately upon solve
        try {
          const raw = localStorage.getItem('nodeword_state_v1');
          const s = raw ? JSON.parse(raw) : {};
          const today = (()=>{const d=new Date();return `${d.getFullYear()}-${d.getMonth()+1}-${d.getDate()}`;})();
          if (s.lastDay !== today) { s.today = 0; s.lastDay = today; }
          s.total = (s.total||0)+1; s.today=(s.today||0)+1; s.best = Math.max(s.best||0, s.today||0);
          s.puzzle = {...(s.puzzle||{}), solved:true};
          localStorage.setItem('nodeword_state_v1', JSON.stringify(s));
          const statTotal = document.getElementById('statTotal'); if (statTotal) statTotal.textContent = String(s.total||0);
          const statToday = document.getElementById('statToday'); if (statToday) statToday.textContent = String(s.today||0);
          const statBest = document.getElementById('statBest'); if (statBest) statBest.textContent = String(s.best||0);
        } catch {}
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
  function readState() {
    try { const raw = localStorage.getItem(STORAGE_KEY); return raw ? JSON.parse(raw) : null; } catch { return null; }
  }
  function canPersist() { return appState && appState.consent === true; }
  function writeState(s) {
    if (!canPersist()) return;
    try { localStorage.setItem(STORAGE_KEY, JSON.stringify(s)); } catch {}
  }
  function todayStamp() {
    const d = new Date();
    return `${d.getFullYear()}-${d.getMonth()+1}-${d.getDate()}`;
  }
  function initState() {
    const today = todayStamp();
    const s = readState();
    if (!s) return { consent: null, total: 0, today: 0, best: 0, lastDay: today, target: 5, puzzle: null };
    if (s.lastDay !== today) { s.today = 0; s.lastDay = today; s.target = 5; s.puzzle = null; }
    return s;
  }
  let appState = initState();
  function updateStatsUI() {
    if (statTotal) statTotal.textContent = String(appState.total || 0);
    if (statToday) statToday.textContent = String(appState.today || 0);
    if (statBest) statBest.textContent = String(appState.best || 0);
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
  // Close stats when clicking backdrop
  statsModal?.addEventListener('click', (e) => { if (e.target === statsModal) closeStatsModal(); });
  // Ensure stats modal starts hidden
  if (statsModal) statsModal.hidden = true;

  // Progressive puzzle size: start at 5, increase to max 12 after each solve (persisted)
  let currentTargetWords = appState.target || 5;
  const maxTargetWords = 12;

  let wordData = null;

  async function generateAndRender() {
    container.textContent = 'Generating…';
    console.log('[Nodeword] Generating puzzle…');
    // Hide next button until solved and clear status message
    if (nextBtn) nextBtn.style.visibility = 'hidden';
    const statusEl = document.getElementById('status');
    if (statusEl) statusEl.textContent = '';
    try {
      if (!wordData) wordData = await fetchWordData();
      const cfg = NODEWORD_CONFIG;
      const target = Math.min(Math.max(5, currentTargetWords), maxTargetWords);
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

      const savedValid = validateSavedPuzzle(appState.puzzle) && appState.puzzle.target === target && isSolvableByCanonical(appState.puzzle.graph);
      let graph;
      if (savedValid) {
        graph = appState.puzzle.graph;
      } else {
        let attempts = 0;
        do {
          if (attempts === 0) console.log('[Nodeword] Generating graph with constraints…', { ...cfg, targetWords: target });
          graph = generatePuzzleGraph(wordData, target, cfg.maxDegree, cfg.maxAttempts);
          attempts++;
        } while (!isSolvableByCanonical(graph) && attempts < 6);
        if (appState.puzzle) { appState.puzzle = null; writeState(appState); }
      }
      console.log('[Nodeword] Graph generated:', {
        words: graph.words.length,
        categories: graph.categories.length,
        edges: graph.edges.length
      });
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
          appState.puzzle = {
            startedAt: appState.puzzle?.startedAt || Date.now(),
            target,
            graph,
            ...payload,
          };
          appState.target = target;
          writeState(appState);
        }
      };
      renderForceGraph(container, aliasGraph, wordData, catEmojis, persist);
      console.log('[Nodeword] Graph rendered to SVG');
    } catch (err) {
      container.textContent = 'Failed to generate puzzle.';
      console.error(err);
    }
  }

  nextBtn?.addEventListener('click', () => {
    console.log('[Nodeword] Next puzzle clicked');
    // Stats increment on solve, and increase target words up to max, then re-generate
    appState.total = (appState.total || 0) + 1;
    appState.today = (appState.today || 0) + 1;
    appState.best = Math.max(appState.best || 0, appState.today || 0);
    writeState(appState);
    updateStatsUI();
    currentTargetWords = Math.min((appState.target || currentTargetWords) + 1, maxTargetWords);
    appState.target = currentTargetWords;
    // Clear current puzzle restore state for next puzzle
    appState.puzzle = null;
    writeState(appState);
    generateAndRender();
  });
  if (nextBtn) nextBtn.style.visibility = 'hidden';
  generateAndRender();
});

