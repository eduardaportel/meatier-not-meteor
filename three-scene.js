export function initThreeScene() {
    // Scene setup
    const scene = new THREE.Scene();
    const camera = new THREE.PerspectiveCamera(75, window.innerWidth / window.innerHeight, 0.1, 1000);
    const renderer = new THREE.WebGLRenderer();
    renderer.setSize(window.innerWidth, window.innerHeight);

    // Add lighting
    const ambientLight = new THREE.AmbientLight(0xffffff, 0.2);
    scene.add(ambientLight);
    const pointLight = new THREE.PointLight(0xffffff, 1.5, 1000);
    scene.add(pointLight);

    // Create Sun
    const sunGeometry = new THREE.SphereGeometry(5, 32, 32);
    const sunMaterial = new THREE.MeshBasicMaterial({ color: 0xfdb813 });
    const sun = new THREE.Mesh(sunGeometry, sunMaterial);
    scene.add(sun);
    
    // Create Earth
    const earthGeometry = new THREE.SphereGeometry(1, 32, 32);
    const earthMaterial = new THREE.MeshStandardMaterial({ color: 0x2e6da4 });
    const earth = new THREE.Mesh(earthGeometry, earthMaterial);
    earth.position.x = 20;
    scene.add(earth);

    // Camera position
    camera.position.z = 40;
    camera.position.y = 10;
    camera.lookAt(scene.position);

    // Animation loop
    const animate = () => {
        requestAnimationFrame(animate);

        // Earth's orbit
        earth.position.x = 20 * Math.cos(Date.now() * 0.0001);
        earth.position.z = -20 * Math.sin(Date.now() * 0.0001);
        earth.rotation.y += 0.005;
        sun.rotation.y += 0.001;

        renderer.render(scene, camera);
    };

    return { domElement: renderer.domElement, animate };
}