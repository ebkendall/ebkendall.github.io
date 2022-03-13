function make2DArray(cols, rows) {
	let arr = new Array(cols);
	for (let i = 0; i < arr.length; i++) {
		arr[i] = new Array(rows);
	}
	return(arr)
}

let grid;
let next;
let cols;
let rows;
let resolution = 10;
let fillCol = [50, 255];

// Switch
let turnOff;
let bottomSame;

// Button selection
let button; 		 // "turnOff" switch
let button_random;
let button_glider;
let button_penta;
let button_gun;
let button_blank;



function setup() {
	createCanvas(600, 450);
	cols = width / resolution;
	rows = height / resolution;

	grid = make2DArray(cols, rows);

	for (let i = 0; i < cols; i++) {
		for (let j = 0; j < rows; j++) {
			grid[i][j] = new Cell(i * resolution, j * resolution);
			grid[i][j].val = floor(random(2));
		}
	}

	next = grid;

	turnOff = 1;
	bottomSame = 0;

	// button characteristics
	button = createButton('START/STOP');
	button.mousePressed(changeB);

	button_glider = createButton('Ex: Glider');
	button_glider.mousePressed(exampleFig1);

	button_penta = createButton('Ex: Pentadecathlon');
	button_penta.mousePressed(exampleFig2);

	button_gun = createButton('Ex: Gosper glider gun');
	button_gun.mousePressed(exampleFig3);

	button_random = createButton('Random');
	button_random.mousePressed(reRandom);

	button_blank = createButton('Blank');
	button_blank.mousePressed(reSet);

}

function draw() {
	background(0);

	for (let i = 0; i < cols; i++) {
		for (let j = 0; j < rows; j++) {
			grid[i][j].display();
		}
	}

	if(turnOff == 0) {
		next = make2DArray(cols, rows);

		for (let i = 0; i < cols; i++) {
			for (let j = 0; j < rows; j++) {

				let state = grid[i][j].val;

				let sum = 0;
				let neighbors = countNeighbors(grid, i , j);

				if (state == 0 && neighbors == 3) {
					next[i][j] = new Cell(i * resolution, j * resolution);
					next[i][j].val = 1;
				} else if (state == 1 && (neighbors < 2 || neighbors > 3)) {
					next[i][j] = new Cell(i * resolution, j * resolution);
					next[i][j].val = 0;
				} else {
					next[i][j] = new Cell(i * resolution, j * resolution);
					next[i][j].val = state;
				}
			}
		}

		if (bottomSame == 1) {
			for(let i = 0; i < cols; i++) {
				next[i][0].val = 0;
			}
		}
		grid = next;
	}
}

function countNeighbors(grid, x, y) {

	let sum = 0;

	for (let i = -1; i < 2; i ++) {
		for (let j = -1; j < 2; j++) {

			let w = (x + cols + i) % cols;
			let z = (y + rows + j) % rows;

			sum += grid[w][z].val;
		}
	}
	sum -= grid[x][y].val;
	return sum;
}

function changeB() {
	if (turnOff == 0) {
		turnOff = 1;
	} else {
		turnOff = 0;
	}
}

function reRandom() {

	turnOff = 1; // Stop the motion and re-randomize things
	bottomSame = 0;

	for (let i = 0; i < cols; i++) {
		for (let j = 0; j < rows; j++) {
			grid[i][j].val = floor(random(2));
		}
	}
}

// Glider
function exampleFig1() {

	turnOff = 1;
	bottomSame = 0;

	// Reset everything to black
	for (let i = 0; i < cols; i++) {
		for (let j = 0; j < rows; j++) {
			grid[i][j].val = 0;
		}
	}

	grid[14][10].val = 1; grid[17][10].val = 1;
	grid[14][12].val = 1; grid[15][13].val = 1;
	grid[16][13].val = 1; grid[17][13].val = 1;
	grid[18][13].val = 1; grid[18][12].val = 1;
	grid[18][11].val = 1;
}

// Penta
function exampleFig2() {

	turnOff = 1;
	bottomSame = 0;

	// Reset everything to black
	for (let i = 0; i < cols; i++) {
		for (let j = 0; j < rows; j++) {
			grid[i][j].val = 0;
		}
	}
	grid[17][10].val = 1; grid[17][11].val = 1; grid[17][13].val = 1;
	grid[17][14].val = 1; grid[17][15].val = 1; grid[17][16].val = 1;
	grid[17][18].val = 1; grid[17][19].val = 1;
	grid[16][12].val = 1; grid[16][17].val = 1;
	grid[18][12].val = 1; grid[18][17].val = 1;

}

// Gun
function exampleFig3() {

	turnOff = 1;
	bottomSame = 1;

	// Reset everything to black
	for (let i = 0; i < cols; i++) {
		for (let j = 0; j < rows; j++) {
			grid[i][j].val = 0;
		}
	}

	// Square
	grid[3][13].val = 1; grid[4][13].val = 1;
	grid[3][14].val = 1; grid[4][14].val = 1;

	// Center
	grid[13][13].val = 1; grid[13][14].val = 1; grid[13][15].val = 1;
	grid[14][12].val = 1; grid[14][16].val = 1;
	grid[15][17].val = 1; grid[16][17].val = 1;
	grid[15][11].val = 1; grid[16][11].val = 1; grid[17][14].val = 1;
	grid[18][12].val = 1; grid[18][16].val = 1;
	grid[19][13].val = 1; grid[19][14].val = 1; grid[19][15].val = 1;
	grid[20][14].val = 1;

	// other V
	grid[23][13].val = 1; grid[23][12].val = 1; grid[23][11].val = 1;
	grid[24][13].val = 1; grid[24][12].val = 1; grid[24][11].val = 1;
	grid[25][14].val = 1; grid[25][10].val = 1;
	grid[27][10].val = 1; grid[27][9].val = 1; grid[27][14].val = 1; grid[27][15].val = 1;

	// Square
	grid[37][11].val = 1; grid[38][11].val = 1;
	grid[37][12].val = 1; grid[38][12].val = 1;

}

// Blank
function reSet() {
	turnOff = 1;
	bottomSame = 0;

	// Reset everything to black
	for (let i = 0; i < cols; i++) {
		for (let j = 0; j < rows; j++) {
			grid[i][j].val = 0;
		}
	}
}

function mousePressed() {
	for (let i = 0; i < cols; i++) {
		for (let j = 0; j < rows; j++) {
			grid[i][j].clicked(mouseX, mouseY);
		}
	}
}

class Cell {
	constructor(x, y) {
		this.x = x;
		this.y = y;
		this.val = 0;
	}

	clicked(x, y) {
		if((x > this.x) && (x < this.x + resolution) &&
		   (y > this.y) && (y < this.y + resolution)) {
			   if(this.val == 1) {
				   this.val = 0;
			   } else {
				   this.val = 1;
			   }
		}

	}

	display() {
		fill(fillCol[this.val]);
		stroke(0);
		rect(this.x, this.y, resolution-1, resolution-1);
	}
}

// I want to add the following:
// (1) Menu for starting screens (only appears during STOP)
// (2) Blank template
// (3) Random template
