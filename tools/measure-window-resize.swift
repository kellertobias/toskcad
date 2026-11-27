// Developer-only macOS window-resize driver (CAD-420). It never ships, links
// or loads into the application. It resizes the main window of a running
// ToskCAD process and prints the epoch times that
// tools/measure-window-resize.mjs uses to cut the response trace.
//
//   swiftc -O tools/measure-window-resize.swift -o "$TMP/resize-driver"
//   "$TMP/resize-driver" PID grip|steps [CYCLES] [DX] [DY] [STEP_MS]
//
// grip:  drags the application's own bottom-right resize grip with
//        synthetic mouse events (it moves the cursor, then restores it).
// steps: sets the window size through Accessibility in small steps, the
//        pattern of window tiling and zoom tools, without moving the cursor.
// The calling terminal needs Accessibility permission for both modes.
import Cocoa

let args = CommandLine.arguments
guard args.count >= 3, let pid = Int32(args[1]), ["grip", "steps"].contains(args[2]) else {
    FileHandle.standardError.write("usage: resize-driver PID grip|steps [CYCLES] [DX] [DY] [STEP_MS]\n".data(using: .utf8)!)
    exit(2)
}
let mode = args[2]
let cycles = args.count > 3 ? Int(args[3])! : 2
let dx = args.count > 4 ? Double(args[4])! : 300
let dy = args.count > 5 ? Double(args[5])! : 150
let stepMs = args.count > 6 ? UInt32(args[6])! : 8
let steps = 60

func attribute(_ element: AXUIElement, _ name: String) -> CFTypeRef? {
    var value: CFTypeRef?
    AXUIElementCopyAttributeValue(element, name as CFString, &value)
    return value
}
let application = AXUIElementCreateApplication(pid)
guard let window = (attribute(application, kAXMainWindowAttribute) ?? attribute(application, kAXFocusedWindowAttribute)) else {
    FileHandle.standardError.write("no window for pid \(pid)\n".data(using: .utf8)!)
    exit(1)
}
let element = window as! AXUIElement
var origin = CGPoint.zero, size = CGSize.zero
AXValueGetValue(attribute(element, kAXPositionAttribute) as! AXValue, .cgPoint, &origin)
AXValueGetValue(attribute(element, kAXSizeAttribute) as! AXValue, .cgSize, &size)
print("window x=\(origin.x) y=\(origin.y) width=\(size.width) height=\(size.height)")

func fraction(_ index: Int, _ outward: Bool) -> Double {
    outward ? Double(index) / Double(steps) : 1 - Double(index) / Double(steps)
}
func stamp(_ label: String) {
    print(String(format: "%@ %.6f", label, Date().timeIntervalSince1970))
}

if mode == "grip" {
    // Inside the 22 pt client grip, clear of AppKit's own resize border.
    let start = CGPoint(x: origin.x + size.width - 11, y: origin.y + size.height - 11)
    let saved = CGEvent(source: nil)!.location
    func post(_ type: CGEventType, _ point: CGPoint) {
        CGEvent(mouseEventSource: nil, mouseType: type, mouseCursorPosition: point, mouseButton: .left)!
            .post(tap: .cghidEventTap)
    }
    post(.mouseMoved, start)
    usleep(150_000)
    stamp("resize-start")
    post(.leftMouseDown, start)
    usleep(30_000)
    for _ in 0..<cycles {
        for outward in [true, false] {
            for index in 1...steps {
                let f = fraction(index, outward)
                post(.leftMouseDragged, CGPoint(x: start.x + dx * f, y: start.y + dy * f))
                usleep(stepMs * 1000)
            }
        }
    }
    post(.leftMouseUp, start)
    stamp("resize-end")
    post(.mouseMoved, saved)
} else {
    stamp("resize-start")
    for _ in 0..<cycles {
        for outward in [true, false] {
            for index in 1...steps {
                let f = fraction(index, outward)
                var next = CGSize(width: size.width + dx * f, height: size.height + dy * f)
                AXUIElementSetAttributeValue(element, kAXSizeAttribute as CFString, AXValueCreate(.cgSize, &next)!)
                usleep(stepMs * 1000)
            }
        }
    }
    stamp("resize-end")
}
