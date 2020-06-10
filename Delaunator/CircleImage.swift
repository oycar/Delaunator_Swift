//
//  CircleImage.swift
//  showMe
//
//  Created by Z Chameleon on 26/5/20.
//  Copyright © 2020 Z Cha. All rights reserved.
//

import SwiftUI

struct CircleImage: View {
  var image: Image
    var body: some View {
      image
        .clipShape(Circle())
        .overlay(
          Circle().stroke(Color.white, lineWidth: 4))
        .shadow(radius:10)
    }
}

struct CircleImage_Previews: PreviewProvider {
    static var previews: some View {
      CircleImage(image: zoneData[0].image)
    }
}
