export function initD3Scene() {
    const width = window.innerWidth;
    const height = window.innerHeight;

    // Create SVG container
    const svg = d3.select("#d3-container")
        .append("svg")
        .attr("width", width)
        .attr("height", height);

    // Add a background representing Earth's surface (simple blue)
    svg.append("rect")
        .attr("width", width)
        .attr("height", height)
        .attr("fill", "#63a0d4");

    // Add text instructions
    svg.append("text")
        .attr("x", width / 2)
        .attr("y", 50)
        .attr("text-anchor", "middle")
        .attr("font-size", "24px")
        .attr("fill", "white")
        .text("Click to simulate an asteroid impact");

    // Click event for impact simulation
    svg.on("click", (event) => {
        const [x, y] = d3.pointer(event);
        
        // Simulate crater
        svg.append("circle")
            .attr("cx", x)
            .attr("cy", y)
            .attr("r", 0)
            .attr("fill", "orange")
            .attr("stroke", "red")
            .attr("stroke-width", 2)
            .transition()
            .duration(500)
            .attr("r", 30 + Math.random() * 40); // Random crater size
    });

    return { svgElement: svg.node() };
}