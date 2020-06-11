# Delaunator_Swift
A port of the Delaunator javascript library of mapbox to the Swift programming language.

This library produces Delaunay triangulations (or the dual Voronoi Mesh) for a set of points defined in an input JSON file; included is a bare bones SwiftUI renderer which will show these structures.


![Example Delaunay Triangulation](./Images/Ukraine_Delaunay.png)

When a triangulation includes near degenrate triangles, i.e. almost collinear, the
points which define the Voronoi diverge will diverge to infinity; making rendering the
Voronoi mesh numerically unstable. But for less extreme triangulations the Voronoi Mesh
will be well behaved.

![Example Delaunay Triangulation](./Images/Ukraine_Voronoi.png)

# Example

The bare minimum usage example would be (just open Delaunator_Sw)
```
// Delaunator_Swift always uses point coordinates of type <Double>
let list = [Point(x:0.0, y:0.0), Point(x:1.0, y:0.0), Point(x:1.0, y:1.0), Point(x:0.0, y:1.0)] // [{x 0, y 0}, {x 1, y 0}, {x 1, y 1}, {x 0, y 1}]

// Calculate the triangulation
var delaunay = Delaunator_Swift(from: list)
delaunay.triangulate()

// The output
print(delaunay.triangles) // [0, 2, 1, 0, 3, 2]\n"
print(delaunay.halfEdges) // "[5, -1, -1, -1, -1, 0]\n"
print(delaunay.hull)      // "[0, 3, 2, 1]\n"
```
```

// Delaunator_Swift always uses point coordinates of type <Double>
let points:Array<Double> = [[168, 180], [168, 178], [168, 179], [168, 181], [168, 183], ...]

// Calculate the triangulation
var delaunay = Delaunator_Swift(from: points)

// Print the triangulation
print(delaunay.triangles)
// [623, 636, 619,  636, 444, 619, ...]]]
```

# Install

If you have Swift5.2 you just need the ```Delaunator.swift``` file for basic
usage.

# Extras

Most of the files are for efficient input, output and a very simple rendering app
for SwiftUI. Data can be input using JSON formatted files; if you want to do this
but don't need the renderer you will need the ```Zone.swift``` file too.



Apart from defining the points to be triangulated programatically data can be
provided using a JSON data structure; which is parsed using this structure

```
// We want to write out a JSON data file
struct StoredZones : Codable {
  let zones: [Zone]

  struct Zone: Codable {
    var name: String
    var points: Array<Array<Double>>?

    // The header
    var header: Header?
  }

  struct Header: Hashable, Codable {
    // Extra stuff
    var action:String?
    var tolerance:Double = Static.Tolerance
    var comment:String?
    var hull: Array<Int>?
    var scale:Double?
  }
}
```

Thus a JSON input file might look like - with one or more zones being defined
each with a distinct point set for triangulation.

```
{
  "zones": [
    {
      "name": "issue13",
      "points": [
        [
          4,
          1
        ],
        [
          3.7974166882130675,
          2.0837249985614585
        ],
        [
          3.2170267516619773,
          3.0210869309396715
        ],
        [
          2.337215067329615,
          3.685489874065187
        ],
        [
          1.276805078389906,
          3.9872025288851036
        ],
        [
          0.17901102978375127,
          3.885476929518457
        ],
        [
          -0.8079039091377689,
          3.3940516818407187
        ],
        [
          -1.550651407188842,
          2.5792964886320684
        ],
        [
          -1.9489192990517052,
          1.5512485534497125
        ],
        [
          -1.9489192990517057,
          0.44875144655029087
        ],
        [
          -1.5506514071888438,
          -0.5792964886320653
        ],
        [
          -0.8079039091377715,
          -1.394051681840717
        ],
        [
          0.17901102978374794,
          -1.8854769295184561
        ],
        [
          1.276805078389902,
          -1.987202528885104
        ],
        [
          2.337215067329611,
          -1.6854898740651891
        ],
        [
          3.217026751661974,
          -1.021086930939675
        ],
        [
          3.7974166882130653,
          -0.08372499856146409
        ]
      ]
    }
  ]
}
```



Work Flo
Yaml -> JSON input
