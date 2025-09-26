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

  for (let attempt = 0; attempt < maxAttempts; attempt++) {
    // Reset graph state each attempt
    wordSet.clear(); categorySet.clear(); edges.clear(); wordDegree.clear(); categoryDegree.clear();
    currentWordToCats.clear(); currentCatToWords.clear();

    const startWord = choice(allWords);
    const startCats = wordToCategories[startWord];
    if (!startCats || startCats.length === 0) continue;
    const firstCategory = choice(startCats);
    if (!addEdge(startWord, firstCategory)) continue;

    // Greedily grow the graph while keeping degrees <= maxDegree and maintaining connectivity
    while (wordSet.size < targetWordCount) {
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
          if (addEdge(w, category)) {
            foundNewWord = true;
            expanded = true;
            break;
          }
        }

        if (foundNewWord) break;
        // rollback edge anchor->category if it didn't lead to expansion
        const key = anchorWord + "||" + category;
        if (edges.has(key)) {
          edges.delete(key);
          wordDegree.set(anchorWord, (wordDegree.get(anchorWord) || 1) - 1);
          categoryDegree.set(category, (categoryDegree.get(category) || 1) - 1);
          // Update adjacency maps
          const wc = currentWordToCats.get(anchorWord);
          if (wc) {
            wc.delete(category);
            if (wc.size === 0) currentWordToCats.delete(anchorWord);
          }
          const cw = currentCatToWords.get(category);
          if (cw) {
            cw.delete(anchorWord);
            if (cw.size === 0) currentCatToWords.delete(category);
          }
          if (!currentCatToWords.has(category)) categorySet.delete(category);
          if (!currentWordToCats.has(anchorWord)) wordSet.delete(anchorWord);
        }
      }

      if (!expanded) break;
      if (!isContiguous()) break;
    }

    if (wordSet.size >= targetWordCount && isContiguous()) {
      // Build final structure
      const words = Array.from(wordSet);
      const categories = Array.from(categorySet);
      const edgesList = [];
      for (const key of edges) {
        const [w, c] = key.split('||');
        edgesList.push({ word: w, category: c });
      }
      return { words, categories, edges: edgesList };
    }
  }

  throw new Error("Failed to generate a contiguous graph within constraints");
}

async function fetchWordData() {
  const res = await fetch('/data/words.json');
  if (!res.ok) throw new Error('Failed to load words.json');
  return res.json();
}

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

function renderForceGraph(container, aliasGraph) {
  container.innerHTML = '';
  const svg = d3.select(container).append('svg').attr('class', 'graph-svg');
  const width = container.clientWidth;
  const height = Math.min(720, Math.max(480, Math.floor(window.innerHeight * 0.7)));
  svg.attr('width', width).attr('height', height);

  const nodes = aliasGraph.nodes.map(n => ({ ...n }));
  const links = aliasGraph.links.map(l => ({ ...l }));

  const link = svg.append('g')
    .attr('stroke-linecap', 'round')
    .selectAll('line')
    .data(links)
    .enter()
    .append('line')
    .attr('class', 'link');

  const node = svg.append('g')
    .selectAll('g')
    .data(nodes, d => d.alias)
    .enter()
    .append('g')
    .attr('class', d => `node ${d.type}`)
    .call(d3.drag()
      .on('start', (event, d) => {
        if (!event.active) simulation.alphaTarget(0.3).restart();
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

  node.append('circle')
    .attr('r', radius);

  node.append('text')
    .text(d => d.alias);

  const simulation = d3.forceSimulation(nodes)
    .force('link', d3.forceLink(links).id(d => d.alias).distance(60).strength(0.9))
    .force('charge', d3.forceManyBody().strength(-220))
    .force('collide', d3.forceCollide().radius(d => radius(d) + 2).strength(0.9))
    .force('center', d3.forceCenter(width / 2, height / 2))
    .force('x', d3.forceX(width / 2).strength(0.05))
    .force('y', d3.forceY(height / 2).strength(0.05));

  simulation.on('tick', () => {
    link
      .attr('x1', d => d.source.x)
      .attr('y1', d => d.source.y)
      .attr('x2', d => d.target.x)
      .attr('y2', d => d.target.y);
    node
      .attr('transform', d => `translate(${d.x},${d.y})`);
  });

  // Resize handling for responsiveness
  const onResize = () => {
    const w = container.clientWidth;
    const h = Math.min(720, Math.max(480, Math.floor(window.innerHeight * 0.7)));
    svg.attr('width', w).attr('height', h);
    simulation.force('center', d3.forceCenter(w / 2, h / 2));
    simulation.alpha(0.3).restart();
  };
  window.addEventListener('resize', onResize, { passive: true });
}

document.addEventListener("DOMContentLoaded", () => {
  const title = document.querySelector("h1");
  if (title) {
    title.animate(
      [ { transform: "translateY(8px)", opacity: 0 }, { transform: "translateY(0)", opacity: 1 } ],
      { duration: 600, easing: "cubic-bezier(.2,.8,.2,1)", fill: "both" }
    );
  }

  const container = document.getElementById('puzzle');
  const nextBtn = document.getElementById('nextBtn');

  let wordData = null;

  async function generateAndRender() {
    container.textContent = 'Generating…';
    try {
      if (!wordData) wordData = await fetchWordData();
      const graph = generatePuzzleGraph(wordData, 12, 4, 6000);
      const aliasGraph = toAliasGraph(graph);
      renderForceGraph(container, aliasGraph);
    } catch (err) {
      container.textContent = 'Failed to generate puzzle.';
      console.error(err);
    }
  }

  nextBtn?.addEventListener('click', generateAndRender);
  generateAndRender();
});

