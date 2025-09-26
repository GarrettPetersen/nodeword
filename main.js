document.addEventListener("DOMContentLoaded", () => {
  const title = document.querySelector("h1");
  if (!title) return;

  // Subtle entrance animation
  title.animate(
    [
      { transform: "translateY(8px)", opacity: 0 },
      { transform: "translateY(0)", opacity: 1 }
    ],
    { duration: 600, easing: "cubic-bezier(.2,.8,.2,1)", fill: "both" }
  );

  // Placeholder: wire up Start button when the game is implemented
  const startButton = document.querySelector("button.primary");
  if (startButton) {
    startButton.addEventListener("click", () => {
      console.log("Nodeword starting... (placeholder)");
    });
  }
});


