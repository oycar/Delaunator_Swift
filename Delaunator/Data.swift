//
//  Data.swift
//  showMe
//
//  Created by Z Chameleon on 26/5/20.
//  Copyright © 2020 Z Cha. All rights reserved.
//


import UIKit
import SwiftUI
import CoreLocation

// Read input data using JSON

let zoneData: [Zone] = readZones(using: load("zoneData.json"))
let triangulationData: [Triangulation] = zoneData.map({$0.triangulation ?? Triangulation()})

func load<T: Decodable>(_ filename: String) -> T {
  let data: Data
  
  guard let file = Bundle.main.url(forResource: filename, withExtension: nil)
    else {
      fatalError("Couldn't find \(filename) in main bundle.")
  }
  
  do {
    data = try Data(contentsOf: file)
  } catch {
    fatalError("Couldn't load \(filename) from main bundle:\n\(error)")
  }
  
  do {
    let decoder = JSONDecoder()
    return try decoder.decode(T.self, from: data)
  } catch {
    fatalError("Couldn't parse \(filename) as \(T.self):\n\(error)")
  }
}

final class ImageStore {
  typealias _ImageDictionary = [String: CGImage]
  fileprivate var images: _ImageDictionary = [:]
  
  fileprivate static var scale = 2
  
  static var shared = ImageStore()
  
  func image(name: String) -> Image {
    let index = _guaranteeImage(name: name)
    
    return Image(images.values[index], scale: CGFloat(ImageStore.scale), label: Text(name))
  }
  
  static func loadImage(name: String) -> CGImage? {
    guard
      let url = Bundle.main.url(forResource: name, withExtension: "jpg"),
      let imageSource = CGImageSourceCreateWithURL(url as NSURL, nil),
      let image = CGImageSourceCreateImageAtIndex(imageSource, 0, nil)
      else {
        return nil
        //fatalError("Couldn't load image \(name).jpg from main bundle.")
    }
    return image
  }
  
  fileprivate func _guaranteeImage(name: String) -> _ImageDictionary.Index {
    if let index = images.index(forKey: name) { return index }
    
    images[name] = ImageStore.loadImage(name: name) ?? ImageStore.loadImage(name: "defaultDelaunay")
    //if let index = images.index(forKey: name) { return index }
    //return images.index(forKey: "ukraine")!
    return images.index(forKey: name)!
  }
}

// Add an image to Zone
extension Zone {
  var image: Image {
    ImageStore.shared.image(name: name)
  }
}


struct Data_Previews: PreviewProvider {
  static var previews: some View {
    /*@START_MENU_TOKEN@*/Text("Hello, World!")/*@END_MENU_TOKEN@*/
  }
}
