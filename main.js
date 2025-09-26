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

function renderGraph(container, graph) {
  const words = new Set(graph.edges.map((e) => e.word));
  const categories = new Set(graph.edges.map((e) => e.category));

  container.innerHTML = '';
  const wrapper = document.createElement('div');
  wrapper.className = 'graph';

  const left = document.createElement('div');
  const right = document.createElement('div');

  const hWords = document.createElement('h3'); hWords.textContent = 'Words';
  const hCats = document.createElement('h3'); hCats.textContent = 'Categories';

  const wordList = document.createElement('div'); wordList.className = 'list';
  for (const w of words) {
    const chip = document.createElement('span'); chip.className = 'chip'; chip.textContent = w;
    wordList.appendChild(chip);
  }

  const catList = document.createElement('div'); catList.className = 'list';
  for (const c of categories) {
    const chip = document.createElement('span'); chip.className = 'chip'; chip.textContent = c;
    catList.appendChild(chip);
  }

  left.appendChild(hWords); left.appendChild(wordList);
  right.appendChild(hCats); right.appendChild(catList);

  wrapper.appendChild(left); wrapper.appendChild(right);

  const edgeTitle = document.createElement('h3'); edgeTitle.textContent = 'Connections';
  const edgeList = document.createElement('div'); edgeList.className = 'list';
  for (const e of graph.edges) {
    const tag = document.createElement('span'); tag.className = 'edge'; tag.textContent = `${e.word} ↔ ${e.category}`;
    edgeList.appendChild(tag);
  }

  container.appendChild(wrapper);
  container.appendChild(edgeTitle);
  container.appendChild(edgeList);
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
      renderGraph(container, graph);
    } catch (err) {
      container.textContent = 'Failed to generate puzzle.';
      console.error(err);
    }
  }

  nextBtn?.addEventListener('click', generateAndRender);
  generateAndRender();
});

