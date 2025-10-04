import { initThreeScene } from './three-scene.js';
import { initD3Scene } from './d3-scene.js';

document.addEventListener('DOMContentLoaded', () => {
    const threeContainer = document.getElementById('three-container');
    const d3Container = document.getElementById('d3-container');
    const toggleButton = document.getElementById('toggle-button');

    let is3D = true;

    // Initialize 3D scene
    const { domElement: threeCanvas, animate: animateThree } = initThreeScene();
    threeContainer.appendChild(threeCanvas);
    animateThree();

    // Initialize 2D scene (it's hidden by default)
    initD3Scene();
    
    // Toggle functionality
    toggleButton.addEventListener('click', () => {
        is3D = !is3D;
        
        if (is3D) {
            threeContainer.style.display = 'block';
            d3Container.style.display = 'none';
            toggleButton.textContent = 'Switch to 2D Impact View';
        } else {
            threeContainer.style.display = 'none';
            d3Container.style.display = 'block';
            toggleButton.textContent = 'Switch to 3D Solar System';
        }
    });

    // Handle window resize
    window.addEventListener('resize', () => {
        // You would add resize logic for both scenes here
    });
});